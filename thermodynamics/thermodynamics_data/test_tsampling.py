from sampling.sampling_dependencies import *
from sampling.optGpSampler_sampling import optGpSampler_sampling
cobra_model = load_json_model("test_tsampling.json");
sampling = optGpSampler_sampling(data_dir_I = "/");
sampling.generate_samples(cobra_model=cobra_model,filename_points = "test_tsampling_points.json",filename_warmup = "test_tsampling_warmup.json",solver_id_I = "optGpSampler",n_points_I = 62,n_steps_I = 5000,n_threads_I = 2);
