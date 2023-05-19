# Small Normal Data:

    Code
      dist_data
    Output
      <MeanNormalData>
        Inherits from: <Data>
        Public:
          alpha: function (x) 
          alpha_subsample: function (x, i) 
          clone: function (deep = FALSE) 
          constant_pi_actual: 4.10321835413825e-11
          constant_pi_observed: 1.13309438425407e-11
          d: 1
          data_x: -3.16154851094265 4.68501868258162 1.47276569197799 3.42 ...
          div_alpha: function (x) 
          div_alpha_subsample: function (x, i) 
          fisher_information: function () 
          get_hessian_bound: function () 
          grad_ll: function (x, i) 
          hessian_bound: 1.09083600348895
          hessian_ll_i: function (i) 
          initialize: function (data_x, mu_true = NULL, sigma_true = NULL) 
          inv_lambda: 0.902912222349122
          lambda: 1.10752737115274
          lap_ll: function (x, i) 
          log_pi_actual: function (x) 
          log_pi_observed: function (x) 
          mu_true: 5
          n: 10
          phi: function (x) 
          phi_control_var: function (x) 
          phi_estimator: function (x, i, j) 
          phi_estimator_bounds: function (maximal_distance) 
          pi_actual: function (x) 
          pi_observed: function (x) 
          precalc_control_var: function (x_hat = NULL) 
          print_debug: function () 
          sigma_est: 3.50230906382132
          sigma_true: 5
          total_grad_ll: function (x) 
          total_lap_ll: function (x) 
          x_hat: 3.95322046100224
          x_hat_C: -0.407625240633715
          x_hat_grad_ll_i: list
          x_hat_lap_ll_i: list
          x_hat_norm_2grad_log: 2.77555756156289e-16
          x_hat_total_grad_ll: -1.38777878078145e-16
          x_hat_total_lap_ll: -0.815250481267431
          x_scale: function (x) 
          x_unscale: function (x) 
        Private:
          grad_ll: function (x, idx) 
          lap_ll: function (x, idx) 

# Large Normal Data:

    Code
      dist_data
    Output
      <MeanNormalData>
        Inherits from: <Data>
        Public:
          alpha: function (x) 
          alpha_subsample: function (x, i) 
          clone: function (deep = FALSE) 
          constant_pi_actual: 1.48838260667195e-41
          constant_pi_observed: 1.99929368622575e-64
          d: 1
          data_x: -3.63230970218853 -2.06299626348368 -2.7054468616044 -2. ...
          div_alpha: function (x) 
          div_alpha_subsample: function (x, i) 
          fisher_information: function () 
          get_hessian_bound: function () 
          grad_ll: function (x, i) 
          hessian_bound: 17.2249489986587
          hessian_ll_i: function (i) 
          initialize: function (data_x, mu_true = NULL, sigma_true = NULL) 
          inv_lambda: 9.66084288275677
          lambda: 0.103510636922256
          lap_ll: function (x, i) 
          log_pi_actual: function (x) 
          log_pi_observed: function (x) 
          mu_true: -2
          n: 100
          phi: function (x) 
          phi_control_var: function (x) 
          phi_estimator: function (x, i, j) 
          phi_estimator_bounds: function (maximal_distance) 
          pi_actual: function (x) 
          pi_observed: function (x) 
          precalc_control_var: function (x_hat = NULL) 
          print_debug: function () 
          sigma_est: 1.03510636922256
          sigma_true: 1
          total_grad_ll: function (x) 
          total_lap_ll: function (x) 
          x_hat: -1.9839822340074
          x_hat_C: -46.665942602656
          x_hat_grad_ll_i: list
          x_hat_lap_ll_i: list
          x_hat_norm_2grad_log: 4.44089209850063e-16
          x_hat_total_grad_ll: -2.22044604925031e-16
          x_hat_total_lap_ll: -93.331885205312
          x_scale: function (x) 
          x_unscale: function (x) 
        Private:
          grad_ll: function (x, idx) 
          lap_ll: function (x, idx) 

# Centered Cauchy Example

    Code
      dist_data
    Output
      <CauchyData>
        Inherits from: <Data>
        Public:
          alpha: function (x) 
          alpha_subsample: function (x, i) 
          cauchy_scale: function (x, mu, gamma) 
          clone: function (deep = FALSE) 
          d: 1
          div_alpha: function (x) 
          div_alpha_subsample: function (x, i) 
          gamma: 1
          grad_ll: function (x, i = NA) 
          hessian_bound: 0.25
          initialize: function (mu, gamma) 
          inv_lambda: 4
          lambda: 0.25
          lap_ll: function (x, i = NA) 
          mu: 0
          n: 1
          phi: function (x) 
          phi_bounds_exact: function (x.l, x.u) 
          phi_control_var: function (x) 
          phi_estimator: function (x, i, j) 
          phi_estimator_bounds: function (maximal_distance) 
          pi: function (x) 
          pi_ll: function (x) 
          precalc_control_var: function (x_hat = NULL) 
          total_grad_ll: function (x) 
          total_lap_ll: function (x) 
          x_hat: 0
          x_hat_C: -1
          x_hat_grad_ll_i: list
          x_hat_lap_ll_i: list
          x_hat_norm_2grad_log: 0
          x_hat_total_grad_ll: 0
          x_hat_total_lap_ll: -2
          x_scale: function (x) 
          x_unscale: function (x) 
        Private:
          grad_ll: function (x, idx) 
          lap_ll: function (x, idx) 

# Deviated Cauchy Example

    Code
      dist_data
    Output
      <CauchyData>
        Inherits from: <Data>
        Public:
          alpha: function (x) 
          alpha_subsample: function (x, i) 
          cauchy_scale: function (x, mu, gamma) 
          clone: function (deep = FALSE) 
          d: 1
          div_alpha: function (x) 
          div_alpha_subsample: function (x, i) 
          gamma: 1
          grad_ll: function (x, i = NA) 
          hessian_bound: 0.25
          initialize: function (mu, gamma) 
          inv_lambda: 4
          lambda: 0.25
          lap_ll: function (x, i = NA) 
          mu: 0
          n: 1
          phi: function (x) 
          phi_bounds_exact: function (x.l, x.u) 
          phi_control_var: function (x) 
          phi_estimator: function (x, i, j) 
          phi_estimator_bounds: function (maximal_distance) 
          pi: function (x) 
          pi_ll: function (x) 
          precalc_control_var: function (x_hat = NULL) 
          total_grad_ll: function (x) 
          total_lap_ll: function (x) 
          x_hat: 0
          x_hat_C: -1
          x_hat_grad_ll_i: list
          x_hat_lap_ll_i: list
          x_hat_norm_2grad_log: 0
          x_hat_total_grad_ll: 0
          x_hat_total_lap_ll: -2
          x_scale: function (x) 
          x_unscale: function (x) 
        Private:
          grad_ll: function (x, idx) 
          lap_ll: function (x, idx) 

