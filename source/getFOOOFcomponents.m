function [full_model,ap_fit,pe_fit] = getFOOOFcomponents(f,aperiodic_params,periodic_params)
% Extract components from parameters in fooof results structure defined in fooof_mat
% f - frequency points
% aperiodic_params - array of aperiodic parameters
% periodic_params - array of periodic parameters
% Return power spectra of full model, aperiodic fit, peak fit in linear scale
try
    full_fit = py.fooof.sim.gen.gen_model(py.numpy.array(f), ...
        py.numpy.array(aperiodic_params),py.numpy.array(periodic_params), ...
        pyargs('return_components',true));
    full_fit = cellfun(@double,cell(full_fit),'UniformOutput',0);
catch
    shape = int32(size(periodic_params));
    periodic_params = periodic_params';
    full_fit = py.fooof.sim.gen.gen_model(py.numpy.array(f), ...
        py.numpy.array(aperiodic_params), ...
        py.numpy.reshape(py.numpy.array(periodic_params(:)),shape), ...
        pyargs('return_components',true));
    full_fit = py.numpy.asarray(full_fit);
    full_fit = cellfun(@double,cell(full_fit.tolist()),'UniformOutput',0);
end
full_model = 10.^full_fit{1}(:);
ap_fit = 10.^full_fit{3}(:);
pe_fit = full_model-ap_fit;
end