#!/usr/bin/env run16-groovy
import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.detector.base.DetectorType
import groovyx.gpars.GParsPool
import org.jlab.clas.physics.LorentzVector
import org.jlab.clas.physics.Vector3
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import org.jlab.groot.data.TDirectory
import java.util.concurrent.ConcurrentHashMap
import java.nio.ByteBuffer
import org.jlab.jroot.ROOTFile
import uconn.utils.pid.Candidate.Level
import uconn.utils.pid.stefan.ElectronCandidate
import uconn.utils.pid.stefan.ElectronCandidate.Cut as ECut
import uconn.utils.pid.stefan.ProtonCandidate
import uconn.utils.pid.stefan.ProtonCandidate.Cut as PCut
import uconn.utils.pid.stefan.PhotonCandidate
import uconn.utils.pid.stefan.PhotonCandidate.Cut as GCut
import my.Sugar
import clasqa.QADB


Sugar.enable()
/////////////////////////////////////////////////////////////////////

def beam = LorentzVector.withPID(11,0,0,10.6041)
def target = LorentzVector.withPID(2212,0,0,0)

def elecuts = [ECut.CC_NPHE, ECut.DC_VERTEX, ECut.EC_OUTER_VS_INNER, ECut.EC_SAMPLING, ECut.PID, ECut.DC_FIDUCIAL_REG1, ECut.DC_FIDUCIAL_REG2, ECut.DC_FIDUCIAL_REG3]
def procuts = [PCut.PID, PCut.FORWARD, PCut.DELTA_VZ]
def phocuts = [GCut.PID, GCut.FORWARD]

def isinb = args[0].contains('inb') || args[0].contains('torus-1')
def ismc = args[0].contains("cache") || args[0].contains("gemc")

def suff = isinb ? 'inb' : 'outb'
if(ismc) suff += '.mc'
else suff += '.qa'

def ff = new ROOTFile("lvl2_eppi0.${suff}.root")
def varlist = 'ex:ey:ez:px:py:pz:g1x:g1y:g1z:g2x:g2y:g2z:esec:psec:g1sec:g2sec:run:status'

if(ismc) varlist += ':ex0:ey0:ez0:px0:py0:pz0'
else varlist += ':ihel:dcx1:dcy1:dcz1:g1v:g1w:g2v:g2w'

def tt = ff.makeNtuple('h22','title',varlist)

GParsPool.withPool 12, {
args.eachParallel{fname->
	println(fname)
	QADB qa = new QADB()

	def reader = new HipoReader()
	reader.open(fname)
	def event = new Event()
	def factory = reader.getSchemaFactory()
	def banks = ['RUN::config','REC::Event','REC::Particle','REC::Calorimeter','REC::Cherenkov','REC::Traj','REC::Scintillator'].collect{new Bank(factory.getSchema(it))}

	while(reader.hasNext()) {
		reader.nextEvent(event)
		banks.each{event.read(it)}

		if(banks.every()) {
			def (runb,evb,partb,ecb,ccb,trajb,scb) = banks

			def run = runb.getInt("run",0)
			def evn = runb.getInt("event",0)

			if(ismc || qa.OkForAsymmetry(run, evn))
			if(true) {

				def elecan = ElectronCandidate.getElectronCandidate(0, partb, ecb, ccb, trajb, isinb)

				def procans = []
				if(elecan.iselectron(*elecuts))
					procans = (0..<partb.getRows()).collect{ProtonCandidate.getProtonCandidate(it, partb, trajb, isinb)}.findAll{it.isproton(*procuts)}

				def gcans = []
				if(procans)
					gcans = (1..<partb.getRows()).collect{PhotonCandidate.getPhotonCandidate(it, partb, ecb, isinb)}.findAll{it.isphoton(*phocuts) && it.getPCALsector()}

				if(gcans.size()>1) {
					int iele = elecan.ipart, esec = elecan.getPCALsector()
					def ele = elecan.getLorentzVector()

					def pi0cans = [gcans, gcans].combinations().findAll{it[1].ipart > it[0].ipart}.findAll{g1can,g2can->
						int ig1 = g1can.ipart, ig2 = g2can.ipart
						int g1sec = g1can.getPCALsector(), g2sec = g2can.getPCALsector()

						def g1 = g1can.getLorentzVector()
						def g2 = g2can.getLorentzVector()

						def dgsec = 10
						if(g1sec && g2sec) {
							dgsec = (g2sec-g1sec).abs()
							if(dgsec>3) dgsec = 6-dgsec
						}

						def fgg = (g1+g2).mass()

						return dgsec==0 && esec!=g1sec && esec!=g2sec && g1.e()>0.6 && g2.e()>0.6 && fgg>0.07 && fgg<0.2
					}

					def finals = []
					if(pi0cans) finals = [procans, pi0cans].combinations().findResults{procan, ggcans->
						def (g1can, g2can) = ggcans

						int ipro = procan.ipart, ig1 = g1can.ipart, ig2 = g2can.ipart
						int psec = procan.getDCsector(), g1sec = g1can.getPCALsector(), g2sec = g2can.getPCALsector()

						def pro = procan.getLorentzVector()
						def g1 = g1can.getLorentzVector()
						def g2 = g2can.getLorentzVector()

						def epx = beam+target-ele-pro
						def gg = g1+g2
						def epggx = epx-gg

						def feggx = (beam+target-ele-gg).mass()
						def fgg = gg.mass()
						def fepx = epx.mass2()
						def fmisse = epggx.e()
						def ftheta = epx.vect().theta(gg.vect())
						def dpx = epggx.px().abs()
						def dpy = epggx.py().abs()
						def dpt = Math.sqrt(dpx*dpx + dpy*dpy)
						def dpz = epggx.pz()

						def dphi = Math.toDegrees(gg.phi()-epx.phi()).abs()
	//          def fpz = epggx.pz()

						def ww = (beam+target-ele).mass()
						def q2 = -(beam-ele).mass2()

						def cuts = [ ww>2, q2>2,
							dpx<0.4, dpy<0.4, dpz>-0.75, dpz<1,
							dphi<5,
							//dpx<0.4, dpy<0.4, dpz>-0.5, dpz<0.75,
							//dphi<4,
							]

						if(cuts.every()) {
							def status = 0

							if(elecan.cut_EC_FIDUCIAL(Level.LOOSE))
								status += 1
							if(elecan.cut_EC_FIDUCIAL(Level.TIGHT))
								status += 2
							if(procan.cut_CHI2PID())
								status += 4
							if(procan.cut_DC_FIDUCIAL_REG1() && procan.cut_DC_FIDUCIAL_REG2() && procan.cut_DC_FIDUCIAL_REG3())
								status += 8
							if(g1can.cut_EC_FIDUCIAL(Level.LOOSEST) && g2can.cut_EC_FIDUCIAL(Level.LOOSEST))
								status += 16
							if(g1can.cut_EC_FIDUCIAL(Level.LOOSE) && g2can.cut_EC_FIDUCIAL(Level.LOOSE))
								status += 32
							if(g1can.cut_BETA() && g2can.cut_BETA())
								status += 64

							def vars = [ele.px(),ele.py(),ele.pz(),
								pro.px(),pro.py(),pro.pz(),
								g1.px(),g1.py(),g1.pz(),
								g2.px(),g2.py(),g2.pz(),
								esec,psec,g1sec,g2sec,
								run,status
							]


							if(ismc) {
								def mcb = new Bank(factory.getSchema("MC::Particle"))
								event.read(mcb)

								def gele = LorentzVector.withPID(11,*['px','py','pz'].collect{mcb.getFloat(it,0)})
								def gpro = LorentzVector.withPID(2212,*['px','py','pz'].collect{mcb.getFloat(it,1)})

								vars += [gele.px(), gele.py(), gele.pz(), gpro.px(), gpro.py(), gpro.pz()]

							} else {
								def ihel = evb.getByte("helicity",0)
								def dcx1=procan.getDC1x(), dcy1=procan.getDC1y(), dcz1=procan.getDC1z()
								vars += [ihel, dcx1,dcy1,dcz1, g1can.pcal_lv,g1can.pcal_lw, g2can.pcal_lv,g2can.pcal_lw]
							}

							return vars
						}

						return null
					}

					if(finals) {
						// def npros = (finals.collect{it[1]} as Set).size()
						// def npi0s = (finals.collect{it[2..3]} as Set).size()

						def fill = finals.max{pp->pp[6..11].sum{it*it}}
						tt.fill(*fill)
					}
				}
			}
		}
	}

	reader.close()
}
}

tt.write()
ff.close()


