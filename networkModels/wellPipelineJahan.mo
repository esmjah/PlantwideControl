within networkModels;
model wellPipelineJahan
  replaceable package components = networkModels.networkComponents;
  Modelica.Blocks.Sources.Constant uw(k = 1);
  Modelica.Blocks.Sources.Constant ur(k = 0.04);
  components.wellJahanshahi w;
  components.annulusJahanshahi a;
  components.pipelineJahanshahi p;
  components.gasManifold gm;
  components.separatorTwoPhase sep(constantPressure = 5010000);
equation
  connect(gm.fout, a.fin);
  connect(a.fout, w.inj);
  connect(w.fout, p.fin);
  connect(p.fout, sep.fin);
  connect(uw.y, w.z1);
  connect(ur.y, p.z);
end wellPipelineJahan;
