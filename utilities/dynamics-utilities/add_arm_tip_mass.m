function params = add_arm_tip_mass(m_add,params)

L = [0 0 -params.lowerArmLength];
l = L - params.elbowCOM;
mtot = m_add + params.elbowMass;
params.elbowCOM = (params.elbowCOM*params.elbowMass +m_add*L)/mtot;

params.elbowMass  = mtot;
Ic = reshape(params.elbowRotationalInertia,3,3);
Ic = Ic+m_add*skew(l)*skew(l)';
params.elbowRotationalInertia = reshape(Ic,1,9);

end