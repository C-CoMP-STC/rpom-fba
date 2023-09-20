import numpy as np
from tqdm import tqdm


def get_interpolator(t, y):
    def Y_t(times): return np.interp(np.atleast_1d(times), t, y)
    return Y_t


def rk45(df_dt, y0, tmin, tmax, dt=0.01, terminate_on_error=True, pbar=True, pbar_desc=None, listeners=None):
    def rk45_step(df_dt, y0, dt):
        k1 = df_dt(y0) * dt
        k2 = df_dt(y0 + 0.5 * k1) * dt
        k3 = df_dt(y0 + 0.5 * k2) * dt
        k4 = df_dt(y0 + k3) * dt

        return y0 + (k1 + 2*k2 + 2*k3 + k4)/6

    t_range = np.arange(tmin, tmax, dt)
    result = np.zeros((t_range.size, y0.size))
    result[0, :] = y0
    listener_data = ([listener(y0) for listener in listeners]
                     if listeners is not None else [])

    t_index = range(1, len(t_range))
    for i in tqdm(t_index, pbar_desc) if pbar else t_index:
        try:
            y = rk45_step(df_dt, result[i-1], dt)
            result[i, :] = y

            # Run listeners
            listener_data += ([listener(y) for listener in listeners]
                              if listeners is not None else [])

        except Exception as e:
            if terminate_on_error:
                return t_range[:i], result[:i, :], listener_data
            raise e

    return t_range, result, listener_data
