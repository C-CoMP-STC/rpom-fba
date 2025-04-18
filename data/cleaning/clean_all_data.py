from parameters.od_calibration import CellDensityRegressor
import data.cleaning.clean_drawdown_data as clean_drawdown_data
import data.cleaning.clean_growth_data as clean_growth_data
import data.cleaning.clean_cue_data as clean_cue_data

def main():
    clean_drawdown_data.main()
    clean_cue_data.main()
    # TODO: growth data requires OD calibration, which requires CUE data - so put OD calib into the order here
    clean_growth_data.main()

if __name__ == "__main__":
    main()
