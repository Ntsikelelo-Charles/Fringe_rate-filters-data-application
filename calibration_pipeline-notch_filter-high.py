# python Combine_data_2h_via_uvd.py 

# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 0h --spw high
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw high --lst 0h 
# python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 100 --lst 0h --spw high


python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 1h --spw high
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 1h --spw high
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 1h --spw high
python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw high --lst 1h 
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 100 --lst 1h --spw high


python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 2h --spw high
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 2h --spw high
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 2h --spw high
python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw high --lst 2h 
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 100 --lst 2h --spw high


python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 3h --spw high
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 3h --spw high
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 3h --spw high
python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw high --lst 3h 
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 100 --lst 3h --spw high

# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode gleam --filter_half_width 0.10e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode diffuse --filter_half_width 0.10e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode data --filter_half_width 0.10e-3 --lst 0h --spw high
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_20_mHz --data_type filtered --spw high --lst 0h 
# python calibration_with_notch_filter.py --filter_name Notch_filter_20_mHz --N 100 --lst 0h --spw high

# python filter_data_and_model.py --filter_type Notch_filter_80_mHz --mode gleam --filter_half_width 0.40e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_80_mHz --mode diffuse --filter_half_width 0.40e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_80_mHz --mode data --filter_half_width 0.40e-3 --lst 0h --spw high
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_80_mHz --data_type filtered --spw high --lst 0h 
# python calibration_with_notch_filter.py --filter_name Notch_filter_80_mHz --N 100 --lst 0h --spw high


# python filter_data_and_model.py --filter_type Notch_filter_100_mHz --mode gleam --filter_half_width 0.50e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_100_mHz --mode diffuse --filter_half_width 0.50e-3 --lst 0h --spw high
# python filter_data_and_model.py --filter_type Notch_filter_100_mHz --mode data --filter_half_width 0.50e-3 --lst 0h --spw high
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_100_mHz --data_type filtered --spw high --lst 0h 
# python calibration_with_notch_filter.py --filter_name Notch_filter_100_mHz --N 100 --lst 0h --spw high







