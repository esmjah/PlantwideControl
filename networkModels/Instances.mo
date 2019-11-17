within networkModels;
package Instances
  record well1
    extends networkModels.networkComponents.gasliftWellParameters(P_r=160e5, k1=0.95, M_G_a=0.021, K_alpha=0.947, K_f=0.901, rho_L=769.98, M_G_r_t=0.021, Z_ca=0.927, Z_ct=0.99, GOR_nom=0, Gamma = 1.7); //1.01803
  end well1;

  record well2
    extends networkModels.networkComponents.gasliftWellParameters(P_r=170e5, k1=0.95, M_G_a=0.021, K_alpha=0.9525, K_f=0.8602, rho_L=785.53, M_G_r_t=0.021, Z_ca=0.933, Z_ct=0.99, GOR_nom=0, Gamma = 1.7);
  end well2;
end Instances;
