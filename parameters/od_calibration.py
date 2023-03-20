import os
import pickle
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import KFold


class CellCountRegressor:
    def __init__(self, crossover=None, a=None, b=None, warn_nonlinear=True):
        self.crossover = np.Infinity if crossover is None else crossover
        self.a = 0 if a is None else a
        self.b = 0 if b is None else b

        self.warn_nonlinear = warn_nonlinear

        self.lm = LinearRegression(fit_intercept=False)

    def _format_data(self, data):
        if not isinstance(data, np.ndarray):
            data = np.array(data)
        if len(data.shape) == 1:
            data = data.reshape(-1, 1)
        return data

    def _train_linear(self, x, y):
        self.lm.fit(x, y)

    def _train_curve(self, x, y):
        # Use log-log linear model as initialization
        # for the curve-fitting optimization
        guess_model = LinearRegression()
        guess_model.fit(np.log(x), np.log(y))
        intercept = guess_model.intercept_[0]
        slope = guess_model.coef_[0, 0]

        # Do curve fitting
        def model_form(od, a, b): return a * od**b
        (self.a, self.b), self.pcov = curve_fit(model_form,
                                                x.flatten(),
                                                y.flatten(),
                                                p0=(np.exp(intercept), slope))

    def train(self, x, y):
        x = self._format_data(x)
        y = self._format_data(y)

        linear_mask = x <= self.crossover
        linear_range_x = self._format_data(x[linear_mask])
        linear_range_y = self._format_data(y[linear_mask])

        # train linear model
        if linear_range_x.size > 0:
            self._train_linear(linear_range_x, linear_range_y)
        else:
            # just fit a flat line...idk if there's something
            # smarter to do here
            self._train_linear([[1]], [[0]])

        # train curve-fit model
        self._train_curve(x, y)

    def _predict_linear(self, od):
        return self.lm.predict(od)

    def _predict_curve(self, od):
        return self.a * od**self.b

    def predict(self, od):
        od = self._format_data(od)
        if self.warn_nonlinear and any(od > self.crossover):
            warnings.warn(
                f"Predicting cell count from OD outside the presumed linear range! (OD > {self.crossover:2e})")

        return np.where(od <= self.crossover, self._predict_linear(od), self._predict_curve(od))


def KFold_cell_count_regressor(od, cells, k=5):
    kf = KFold(n_splits=k)

    odmax_domain = np.linspace(0, np.max(od), 100)
    mse_scan = np.zeros_like(odmax_domain, dtype="float")

    min_mse = np.Infinity
    model = CellCountRegressor()
    best_regressor = model

    for i, odmax in enumerate(odmax_domain):
        cv_err = 0
        for train, test in kf.split(od):
            x_in = od[train]
            y_in = cells[train]
            x_out = od[test]
            y_out = cells[test]

            model = CellCountRegressor(odmax, warn_nonlinear=False)
            model.train(x_in, y_in)
            pred = model.predict(x_out)

            try:
                cv_err += mean_squared_error(np.log(y_out + 1),
                                             np.log(pred + 1)) / k
            except:
                print(x_in)
                print(pred)

        if cv_err < min_mse:
            min_mse = cv_err
            best_regressor = model

        mse_scan[i] = cv_err

    # retrain regressor using the crossover identified, on all the data
    best_regressor = CellCountRegressor(best_regressor.crossover)
    best_regressor.train(od, cells)

    return best_regressor, odmax_domain, mse_scan


def get_od_to_cell_count_calibration():
    SAVED_REGRESSOR = "parameters/conversions/od_to_cell_count_params.pickle"
    try:
        with open(SAVED_REGRESSOR, "rb") as f:
            return pickle.load(f)
    except FileNotFoundError:
        main()
        with open(SAVED_REGRESSOR, "rb") as f:
            return pickle.load(f)


def main():
    OUT = "parameters/conversions/od_to_cell_count.pickle"
    DATA = "data/ODcounts_calibration.xlsx"
    raw_data = pd.read_excel(DATA)

    # Clean data
    raw_data['OD600 (1-cm)'] /= raw_data["Scale_OD"]
    raw_data["Mean_cells"] = raw_data["Cells per ml, cyto"]
    clean_data = raw_data.sort_values("OD600 (1-cm)")

    od = clean_data['OD600 (1-cm)'].values.reshape(-1, 1)
    cells = clean_data["Mean_cells"].values.reshape(-1, 1)


    best_reg, _, _ = KFold_cell_count_regressor(od, cells, k=4)

    print("Fit regressor with parameters:")
    print(f"\ta = {best_reg.a}")
    print(f"\tb = {best_reg.b}")
    print(f"\tcrossover = {best_reg.crossover}")
    print(f"\tm = {best_reg.lm.coef_[0, 0]}")

    print(f"saving to {OUT}...")

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "wb") as f:
        pickle.dump(best_reg, f)


if __name__ == "__main__":
    main()
