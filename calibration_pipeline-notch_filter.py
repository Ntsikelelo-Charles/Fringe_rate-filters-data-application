# python filter_data_and_model.py --filter_type Notch_filter_40_mHz_interpolated --mode data --filter_half_width 0.20e-3 --lst 0h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 0h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 0h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw low --lst 0h 
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 0h --spw low

# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 1h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 1h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw low --lst 1h 
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz_interpolated --mode data --filter_half_width 0.20e-3 --lst 1h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 1h --spw low

# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 2h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 2h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw low --lst 2h 
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz_interpolated --mode data --filter_half_width 0.20e-3 --lst 2h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 2h --spw low

# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode gleam --filter_half_width 0.20e-3 --lst 3h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 3h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_40_mHz --data_type filtered --spw low --lst 3h 
# python filter_data_and_model.py --filter_type Notch_filter_40_mHz_interpolated --mode data --filter_half_width 0.20e-3 --lst 3h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 3h --spw low

python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 0h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 1h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 2h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 3h --spw low --mode gleam

python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 0h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 1h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 2h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 3h --spw low --mode full

# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode gleam --filter_half_width 0.10e-3 --lst 0h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode diffuse --filter_half_width 0.10e-3 --lst 0h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_20_mHz --data_type filtered --spw low --lst 0h 
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz_interpolated --mode data --filter_half_width 0.10e-3 --lst 0h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 0h --spw low

# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode gleam --filter_half_width 0.10e-3 --lst 1h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode diffuse --filter_half_width 0.10e-3 --lst 1h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_20_mHz --data_type filtered --spw low --lst 1h 
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz_interpolated --mode data --filter_half_width 0.10e-3 --lst 1h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 1h --spw low


# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode gleam --filter_half_width 0.10e-3 --lst 2h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode diffuse --filter_half_width 0.10e-3 --lst 2h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_20_mHz --data_type filtered --spw low --lst 2h 
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz_interpolated --mode data --filter_half_width 0.10e-3 --lst 2h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 2h --spw low


# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode gleam --filter_half_width 0.10e-3 --lst 3h --spw low
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz --mode diffuse --filter_half_width 0.10e-3 --lst 3h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Notch_filter_20_mHz --data_type filtered --spw low --lst 3h 
# python filter_data_and_model.py --filter_type Notch_filter_20_mHz_interpolated --mode data --filter_half_width 0.10e-3 --lst 3h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 3h --spw low


python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 0h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 1h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 2h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 3h --spw low --mode gleam

python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 0h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 1h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 2h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_20_mHz_interpolated --N 250 --lst 3h --spw low --mode full





