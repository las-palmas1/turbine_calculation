from scipy.integrate import quad


class Calculation:
    def __init__(self, sigma_t11, sigma_r11, omega, theta, theta_arr, r_arr, E_arr, h_arr, rho_new_arr, mu):
        self.omega = omega
        self.sigma_r11 = sigma_r11
        self.theta = theta
        self.theta_arr = theta_arr
        self.r_arr = r_arr
        self.E_arr = E_arr
        self.h_arr = h_arr
        self.mu = mu
        self.rho_new_arr = rho_new_arr
        self.sigma_t_arr = list()
        self.sigma_t_arr.append(list())
        self.sigma_t_arr[0].append(sigma_t11)
        self.sigma_r_arr = list()
        self.sigma_r_arr.append(list())
        self.sigma_r_arr[0].append(sigma_r11)

        for i in range(len(self.r_arr) - 1):

            print(i)

            if i != len(self.r_arr) - 1:
                self.sigma_t_arr.append([])
                self.sigma_r_arr.append([])

            self.sigma_t_1 = self.sigma_t_arr[i][0]
            self.sigma_r_1 = self.sigma_r_arr[i][0]
            self.D_1 = self.sigma_t_1 - self.sigma_r_1
            self.S_1 = self.sigma_t_1 + self.sigma_r_1
            self.fun1_S_2 = 0.5 * (1 + self.mu) * self.rho_new_arr[i] * self.omega ** 2 * (self.r_arr[i + 1] ** 2 -
                                                                                           self.r_arr[i] ** 2)
            self.fun2_S_2 = self.E_arr[i] * (self.theta_arr[i + 1] - self.theta_arr[i])
            self.S_2 = self.S_1 - self.fun1_S_2 - self.fun2_S_2
            self.fun1_D_2 = self.D_1 * (self.r_arr[i] ** 2) / (self.r_arr[i + 1] ** 2)

            def int_D_2(r):
                return self.theta(r) * r

            self.fun2_D_2 = 2 * self.E_arr[i] * quad(int_D_2, self.r_arr[i],
                                                     self.r_arr[i + 1])[0] / (self.r_arr[i + 1] ** 2)
            self.fun3_D_2 = 0.25 * (1 - self.mu) * self.rho_new_arr[i] * self.omega ** 2 * (self.r_arr[i + 1] ** 2 -
                                                                        self.r_arr[i] ** 4 / (self.r_arr[i + 1] ** 2))
            self.fun4_D_2 = self.E_arr[i] * (self.theta_arr[i + 1] -
                                             self.theta_arr[i] * self.r_arr[i] ** 2 / self.r_arr[i + 1] ** 2)
            self.D_2 = self.fun1_D_2 + self.fun2_D_2 + self.fun3_D_2 - self.fun4_D_2
            self.sigma_t_2 = 0.5 * (self.D_2 + self.S_2)
            self.sigma_r_2 = 0.5 * (self.S_2 - self.D_2)
            self.sigma_r_arr[i].append(self.sigma_r_2)
            self.sigma_t_arr[i].append(self.sigma_t_2)
            self.sigma_r_3 = self.sigma_r_2 * self.h_arr[i] / self.h_arr[i + 1]
            self.fun1_s_3 = self.mu * self.sigma_r_2 * self.h_arr[i] / self.h_arr[i + 1]
            self.fun2_s_3 = self.E_arr[i + 1] * (self.sigma_t_2 - self.mu * self.sigma_r_2) / self.E_arr[i]
            self.sigma_t_3 = self.fun1_s_3 + self.fun2_s_3

            if i != len(self.r_arr) - 1:
                self.sigma_t_arr[i + 1].append(self.sigma_t_3)
                self.sigma_r_arr[i + 1].append(self.sigma_r_3)
