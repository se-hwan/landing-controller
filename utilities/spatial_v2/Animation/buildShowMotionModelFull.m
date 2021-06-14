function model = buildShowMotionModelFull(params)

switch params.model
    case 'quad3D'
        model = buildShowMotionModelMC3DFull(params);
    case 'humanoid3D'
        model = buildShowMotionModelH3DFull(params);
end