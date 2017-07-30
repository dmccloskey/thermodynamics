#VSCode remote debugging
######
import ptvsd
ptvsd.enable_attach("my_secret", address = ('0.0.0.0', 3000))
#enable the below line of code only if you want the application to wait until the debugger has attached to it
ptvsd.wait_for_attach()
######

from smartPeak.__main__ import __main__
m = __main__()

# # Test openSWATH_py
# m.run_openSWATH_py(
#     filename_filenames='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_SWATH_filenames.csv',
#     filename_params='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_MRMFeatureFinderScoring_params.csv',
#     delimiter=','
#     )

# # Test openSWATH_cmd
# # filename='/home/user/openMS_MRMworkflow/openSWATH_cmd_params_QC1.csv'
# filename='/home/user/openMS_MRMworkflow/openSWATH_cmd_params_QC1_FeatureXML2TSV.csv'
# m.run_openSWATH_cmd(filename)

# # Test MRMTransitionGroupPicker_py
# m.run_MRMTransitionGroupPicker_py(
#     filename_filenames='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_SWATH_filenames.csv',
#     filename_params='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_MRMFeatureFinderScoring_params.csv',
#     delimiter=','
#     )

# # Test PeakPickerMRM_py
# m.run_PeakPickerMRM_py(
#     filename_filenames='/home/user/openMS_MRMworkflow/BloodProject01_SWATH_filenames.csv',
#     filename_params='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_PeakPickerMRM_params.csv',
#     delimiter=','
#     )

# Test file conversions
# m.convert_MQQMethod2Feature(
#     filename_I='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_qmethod.csv',
#     filename_O='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_SWATH_feature.csv'
#     )

# # Test ReferenceData
# m.run_get_referenceData(
#     experiment_ids_I = ['BloodProject01'],
#     sample_names_I = ['150601_0_BloodProject01_PLT_QC_Broth-1'],
#     acquisition_methods_I = ['140718_McCloskey2013'],
#     used__I = True,
#     settings_filename_I = '/home/user/openMS_MRMworkflow/settings_metabolomics.ini',
#     data_filename_O = '/home/user/openMS_MRMworkflow/BloodProject01/150601_0_BloodProject01_PLT_QC_Broth-1_referenceData.csv')

# Test openSWATH validation
m.run_validate_openSWATH(
    filename_filenames='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_SWATH_filenames.csv',
    filename_params='/home/user/openMS_MRMworkflow/BloodProject01/BloodProject01_validation_params.csv',
    delimiter=','
    )