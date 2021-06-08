function model = buildShowMotionModel(params, model, ground_height)

switch params.model
    case 'quad3D'
        model = buildShowMotionModelMC3D(params, model, ground_height);
    case 'humanoid3D'
        model = buildShowMotionModelHumanoid3D(params, model, ground_height);        
end