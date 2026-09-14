python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 0h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 1h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 2h --spw low --mode gleam
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 3h --spw low --mode gleam

python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 0h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 1h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 2h --spw low --mode full
python FLag_data_post_calibration.py --filter_name Notch_filter_40_mHz_interpolated --N 250 --lst 3h --spw low --mode full

# python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode gleam --filter_half_width 0.20e-3 --lst 0h --spw low
# python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode diffuse --filter_half_width 0.20e-3 --lst 0h --spw low
# python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode data --filter_half_width 0.20e-3 --lst 0h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Main_lobe_baseline_dependent --data_type filtered --spw low --lst 0h
python calibration_with_notch_filter.py --filter_name Main_lobe_baseline_dependent --N 250 --lst 0h --spw low


python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode gleam --filter_half_width 0.20e-3 --lst 1h --spw low
python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode diffuse --filter_half_width 0.20e-3 --lst 1h --spw low
python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode data --filter_half_width 0.20e-3 --lst 1h --spw low
python expand_model_vis_by_redundancy.py --filter_name Main_lobe_baseline_dependent --data_type filtered --spw low --lst 1h
python calibration_with_notch_filter.py --filter_name Main_lobe_baseline_dependent --N 250 --lst 1h --spw low


python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode gleam --filter_half_width 0.20e-3 --lst 2h --spw low
python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode diffuse --filter_half_width 0.20e-3 --lst 2h --spw low
python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode data --filter_half_width 0.20e-3 --lst 2h --spw low
python expand_model_vis_by_redundancy.py --filter_name Main_lobe_baseline_dependent --data_type filtered --spw low --lst 2h
python calibration_with_notch_filter.py --filter_name Main_lobe_baseline_dependent --N 250 --lst 2h --spw low


python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode gleam --filter_half_width 0.20e-3 --lst 3h --spw low
python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode diffuse --filter_half_width 0.20e-3 --lst 3h --spw low
python filter_data_and_model_main_lobe.py --filter_type Main_lobe_baseline_dependent --mode data --filter_half_width 0.20e-3 --lst 3h --spw low
python expand_model_vis_by_redundancy.py --filter_name Main_lobe_baseline_dependent --data_type filtered --spw low --lst 3h
python calibration_with_notch_filter.py --filter_name Main_lobe_baseline_dependent --N 250 --lst 3h --spw low





