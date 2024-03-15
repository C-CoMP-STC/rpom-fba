import numpy as np
from tqdm import tqdm


def get_interpolator(t, y):
    def Y_t(times): return np.interp(np.atleast_1d(times), t, y)
    return Y_t


def runge_kutta(df_dt,
         y0,
         tmin,
         tmax,
         dt=0.01,
         terminate_on_error=True,
         pbar=True,
         pbar_desc=None,
         listeners=None):
    def rk_step(df_dt, y0, dt, t=None):
        k1 = df_dt(y0, t=t) * dt
        k2 = df_dt(y0 + 0.5 * k1, t=t) * dt
        k3 = df_dt(y0 + 0.5 * k2, t=t) * dt
        k4 = df_dt(y0 + k3, t=t) * dt

        return y0 + (k1 + 2*k2 + 2*k3 + k4)/6

    t_range = np.arange(tmin, tmax, dt)

    result = np.zeros((t_range.size, y0.size))
    result[0, :] = y0
    listener_data = ([[listener(y0, tmin)] for listener in listeners]
                     if listeners is not None else [])

    t_index = range(1, len(t_range))
    for i in tqdm(t_index, pbar_desc) if pbar else t_index:
        try:
            y = rk_step(df_dt, result[i-1], dt, t=t_range[i])
            result[i, :] = y

            # Run listeners
            if listeners is not None:
                for l, listener in enumerate(listeners):
                    listener_data[l].append(listener(y, t=t_range[i]))
            # listener_data += ([listener(y, t=t_range[i]) for listener in listeners]
            #                   if listeners is not None else [])

        except Exception as e:
            if terminate_on_error:
                return t_range[:i], result[:i, :], listener_data
            raise e

    return t_range, result, listener_data
