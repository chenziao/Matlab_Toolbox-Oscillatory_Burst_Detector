function PSD_fit = getBackgroundFit(rdat,f)
% Calculate fit background PSD from results structrue.
% results - results structrue
% f - frequency points (use results.f if not specified)
if nargin<2
    f = rdat.f;
end
if isfield(rdat,'aperiodic_params')
    % use fit by FOOOF
    aperiodic_params = py.numpy.array(rdat.aperiodic_params);
    PSD_fit = py.fooof.sim.gen.gen_aperiodic(py.numpy.array(f),aperiodic_params);
    try
        PSD_fit = 10.^double(PSD_fit);
    catch
        PSD_fit = 10.^double(PSD_fit.tolist());
    end
    PSD_fit = PSD_fit(:);
else
    % use linear fit
    PSD_fit = rdat.fit_a*f.^rdat.fit_b;
end
end