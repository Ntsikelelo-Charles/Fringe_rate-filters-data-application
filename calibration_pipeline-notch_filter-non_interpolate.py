
python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 0h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 250 --lst 0h --spw low

python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 1h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 250 --lst 1h --spw low

python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 2h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 250 --lst 2h --spw low

python filter_data_and_model.py --filter_type Notch_filter_40_mHz --mode data --filter_half_width 0.20e-3 --lst 3h --spw low
python calibration_with_notch_filter.py --filter_name Notch_filter_40_mHz --N 250 --lst 3h --spw low











