#!/usr/bin/env run-groovy
import org.jlab.jnp.hipo4.io.HipoReader
import org.jlab.jnp.hipo4.io.HipoWriter
import org.jlab.jnp.hipo4.data.Bank
import org.jlab.jnp.hipo4.data.Event
import org.jlab.jnp.hipo4.data.SchemaFactory
import org.jlab.clas.physics.LorentzVector
import uconn.utils.pid.stefan.ElectronCandidate


import my.Sugar

Sugar.enable()

def beam = LorentzVector.withPID(11,0,0,10.6041)
def targ = LorentzVector.withPID(2212,0,0,0)

def fname = args[0]
def reader = new HipoReader()
reader.open(fname)
def event = new Event()
def factory = reader.getSchemaFactory()
def banknames = ['RUN::config', 'REC::Particle','REC::Calorimeter','REC::Cherenkov','REC::Traj']
def banks = banknames.collect{new Bank(factory.getSchema(it))}

def writer = new HipoWriter(factory)
writer.open('lvl1_eppi0.'+fname.split("/")[-1])

while(reader.hasNext()) {
  reader.nextEvent(event)
  banks.each{event.read(it)}

  if(banks.every()) {
    def (runb,partb,calb,ccb,trajb) = banks

    def pids = (0..<partb.getRows()).findAll{partb.getShort('status',it).div(2000).toInteger().abs()==1}.collect{[partb.getInt('pid',it), it]}.groupBy{it[0]}
    //println(pids)
    if(pids[11] && pids[2212] && pids[22])
    if(pids[22].size()>1) {
      def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,0)})
      def q2 = (beam-ele).mass2();

      def igs = pids[22].collect{it[1]}
      def ispi0 = [igs,igs].combinations().findAll{ig1,ig2 -> ig2>ig1}
        .any{ig1,ig2->
          def g1 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig1)})
          def g2 = LorentzVector.withPID(22,*['px','py','pz'].collect{partb.getFloat(it,ig2)})
          def pi0 = g1+g2;

          return pi0.mass()>0.07 && pi0.mass()<0.2 && q2<-1
        }
      if(ispi0) writer.addEvent(event)
    }
  }
}

reader.close()
writer.close()

