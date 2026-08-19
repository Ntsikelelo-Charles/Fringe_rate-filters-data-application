# python Combine_data_2h_via_uvd.py 
# python Flag_data_for_RFI.py --N 250 --lst 0h 
# python Flag_data_for_RFI.py --N 250 --lst 1h 
# python Flag_data_for_RFI.py --N 250 --lst 2h
# python Flag_data_for_RFI.py --N 250 --lst 3h 

# python expand_model_vis_by_redundancy.py --data_type no_filter --spw low --lst 0h
# python calibration_no_filter.py --N 250 --lst 0h --spw low


# python expand_model_vis_by_redundancy.py --data_type no_filter --spw low --lst 1h
# python calibration_no_filter.py --N 250 --lst 1h --spw low

# python expand_model_vis_by_redundancy.py --filter_name Main_lobe_20_mHz --data_type no_filter --spw low --lst 2h
# python calibration_no_filter.py --N 250 --lst 2h --spw low

python expand_model_vis_by_redundancy.py --data_type no_filter --spw low --lst 3h
python calibration_no_filter.py --N 250 --lst 3h --spw low

python expand_model_vis_by_redundancy.py  --data_type no_filter --spw high --lst 0h
python calibration_no_filter.py --N 250 --lst 0h --spw high

python expand_model_vis_by_redundancy.py --data_type no_filter --spw high --lst 1h
python calibration_no_filter.py --N 250 --lst 1h --spw high

python expand_model_vis_by_redundancy.py --data_type no_filter --spw high --lst 2h
python calibration_no_filter.py --N 250 --lst 2h --spw high

python expand_model_vis_by_redundancy.py --data_type no_filter --spw high --lst 3h
python calibration_no_filter.py --N 250 --lst 3h --spw high

