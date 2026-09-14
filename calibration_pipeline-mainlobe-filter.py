
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode data --filter_half_width 0.20e-3 --lst 0h --spw low
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode gleam --filter_half_width 0.20e-3 --lst 0h --spw low
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 0h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Main_lobe_20_mHz --data_type filtered --spw low --lst 0h
# python calibration_with_notch_filter.py --filter_name Main_lobe_20_mHz --N 250 --lst 0h --spw low

# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode gleam --filter_half_width 0.20e-3 --lst 1h --spw low
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 1h --spw low
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode data --filter_half_width 0.20e-3 --lst 1h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Main_lobe_20_mHz --data_type filtered --spw low --lst 1h
# python calibration_with_notch_filter.py --filter_name Main_lobe_20_mHz --N 250 --lst 1h --spw low


# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode gleam --filter_half_width 0.20e-3 --lst 2h --spw low
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 2h --spw low
# python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode data --filter_half_width 0.20e-3 --lst 2h --spw low
# python expand_model_vis_by_redundancy.py --filter_name Main_lobe_20_mHz --data_type filtered --spw low --lst 2h
# python calibration_with_notch_filter.py --filter_name Main_lobe_20_mHz --N 250 --lst 2h --spw low

python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode gleam --filter_half_width 0.20e-3 --lst 3h --spw low
python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode diffuse --filter_half_width 0.20e-3 --lst 3h --spw low
python filter_data_and_model.py --filter_type Main_lobe_20_mHz --mode data --filter_half_width 0.20e-3 --lst 3h --spw low
python expand_model_vis_by_redundancy.py --filter_name Main_lobe_20_mHz --data_type filtered --spw low --lst 3h
python calibration_with_notch_filter.py --filter_name Main_lobe_20_mHz --N 250 --lst 3h --spw low





