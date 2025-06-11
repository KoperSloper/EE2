# Noise-Robust Realized Variance Project

This project studies how different realized variance estimators perform in the presence of market microstructure noise. It tests the robustness of the ZMA estimator—originally proposed by Zhang, Mykland & Aït‐Sahalia (2005)—under realistic noise settings, and compares it to simpler alternatives like Naïve and 5-minute sampling.

## Files
- **heston.ipynb**  
  Simulates Heston paths, adds microstructure noise (i.i.d., heteroskedastic, AR(1)), and compares Naïve, Five-Minute, Subsampled, and ZMA estimators.

- **Comparisons Naive vs ZMA likelihood.R**  
  Loads the daily realized variances, fits Realized GARCH models in R, and compares out-of-sample log-likelihoods for Naïve vs. ZMA (with Vuong tests).
