use rand::Rng;
use std::f64::consts::E;
use pyo3::prelude::*;

/// Calculates the price of an American call option using the Monte Carlo simulation method.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run.
/// * `num_steps` - The number of steps in each simulation.
///
/// # Returns
///
/// The price of the American call option.
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::monte_carlo_american_call;
///
/// let spot_price = 100.0;
/// let strike_price = 110.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let price = monte_carlo_american_call(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("American call option price: {:.2}", price);
/// ```
#[pyfunction]
fn monte_carlo_american_call(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let dt = time_to_maturity / num_steps as f64;
    let discount_factor = E.powf(-risk_free_rate * time_to_maturity);

    let mut rng = rand::thread_rng();
    let mut payoffs = Vec::with_capacity(num_simulations);

    for _ in 0..num_simulations {
        let mut price = spot_price;
        let mut max_payoff: f64 = 0.0;

        for _ in 0..num_steps {
            let rand_num = rng.gen_range(-1.0..1.0);
            price *= E.powf(
                (risk_free_rate - 0.5 * volatility.powi(2)) * dt
                    + volatility * rand_num * dt.sqrt(),
            );
            max_payoff = max_payoff.max((price - strike_price).max(0.0));
        }

        payoffs.push(max_payoff);
    }

    let sum_payoffs: f64 = payoffs.iter().sum();
    let avg_payoff = sum_payoffs / num_simulations as f64;

    Ok(avg_payoff * discount_factor)
}

/// Calculates the price of an American put option using the Monte Carlo simulation method.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run.
/// * `num_steps` - The number of steps in each simulation.
///
/// # Returns
///
/// The price of the American put option.
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::monte_carlo_american_put;
///
/// let spot_price = 100.0;
/// let strike_price = 90.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let price = monte_carlo_american_put(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("American put option price: {:.2}", price);
/// ```
#[pyfunction]
fn monte_carlo_american_put(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let dt = time_to_maturity / num_steps as f64;
    let discount_factor = E.powf(-risk_free_rate * time_to_maturity);

    let mut rng = rand::thread_rng();
    let mut payoffs = Vec::with_capacity(num_simulations);

    for _ in 0..num_simulations {
        let mut price = spot_price;
        let mut max_payoff: f64 = 0.0;

        for _ in 0..num_steps {
            let rand_num = rng.gen_range(-1.0..1.0);
            price *= E.powf(
                (risk_free_rate - 0.5 * volatility.powi(2)) * dt
                    + volatility * rand_num * dt.sqrt(),
            );
            max_payoff = max_payoff.max((strike_price - price).max(0.0));
        }

        payoffs.push(max_payoff);
    }

    let sum_payoffs: f64 = payoffs.iter().sum();
    let avg_payoff = sum_payoffs / num_simulations as f64;

    Ok(avg_payoff * discount_factor)
}

/// Calculates the delta of an American call option using the finite difference method.
///
/// Delta represents the rate of change of the option price with respect to the change in the underlying asset price.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run for the Monte Carlo method.
/// * `num_steps` - The number of steps in each simulation for the Monte Carlo method.
///
/// # Returns
///
/// The delta of the American call option.
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::{calculate_delta, monte_carlo_american_call};
///
/// let spot_price = 100.0;
/// let strike_price = 110.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let delta = calculate_delta(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("Delta of the American call option: {:.4}", delta);
/// ```
#[pyfunction]
fn calculate_delta(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let h = 0.01;
    let price_up = monte_carlo_american_call(
        spot_price * (1.0 + h),
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    let price_down = monte_carlo_american_call(
        spot_price * (1.0 - h),
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    Ok((price_up - price_down) / (2.0 * h * spot_price))
}

/// Calculates the gamma of an American call option using the finite difference method.
///
/// Gamma represents the rate of change of the option's delta with respect to the change in the underlying asset price.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run for the Monte Carlo method.
/// * `num_steps` - The number of steps in each simulation for the Monte Carlo method.
///
/// # Returns
///
/// The gamma of the American call option.
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::{calculate_gamma, monte_carlo_american_call};
///
/// let spot_price = 100.0;
/// let strike_price = 110.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let gamma = calculate_gamma(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("Gamma of the American call option: {:.4}", gamma);
/// ```
#[pyfunction]
fn calculate_gamma(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let h = 0.01;
    let price_up = monte_carlo_american_call(
        spot_price * (1.0 + h),
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    let price = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    let price_down = monte_carlo_american_call(
        spot_price * (1.0 - h),
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    Ok((price_up - 2.0 * price + price_down) / (h * spot_price).powi(2))
}

/// Calculates the vega of an American call option using the finite difference method.
///
/// Vega represents the sensitivity of the option price to changes in the volatility of the underlying asset.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run for the Monte Carlo method.
/// * `num_steps` - The number of steps in each simulation for the Monte Carlo method.
///
/// # Returns
///
/// The vega of the American call option.
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::{calculate_vega, monte_carlo_american_call};
///
/// let spot_price = 100.0;
/// let strike_price = 110.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let vega = calculate_vega(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("Vega of the American call option: {:.4}", vega);
/// ```
#[pyfunction]
fn calculate_vega(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let h = 0.01;
    let price_up = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate,
        volatility * (1.0 + h),
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    let price_down = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate,
        volatility * (1.0 - h),
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    Ok((price_up - price_down) / (2.0 * h * volatility))
}

/// Calculates the theta of an American call option using the finite difference method.
///
/// Theta represents the rate of change of the option price with respect to the passage of time.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run for the Monte Carlo method.
/// * `num_steps` - The number of steps in each simulation for the Monte Carlo method.
///
/// # Returns
///
/// The theta of the American call option (per day).
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::{calculate_theta, monte_carlo_american_call};
///
/// let spot_price = 100.0;
/// let strike_price = 110.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let theta = calculate_theta(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("Theta of the American call option (per day): {:.4}", theta);
/// ```
#[pyfunction]
fn calculate_theta(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let h = 0.01;
    let price_up = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity * (1.0 + h),
        num_simulations,
        num_steps,
    )?;
    let price_down = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate,
        volatility,
        time_to_maturity * (1.0 - h),
        num_simulations,
        num_steps,
    )?;
    Ok(-(price_up - price_down) / (2.0 * h * time_to_maturity) / 365.0)
}

/// Calculates the rho of an American call option using the finite difference method.
///
/// Rho represents the sensitivity of the option price to changes in the risk-free interest rate.
///
/// # Arguments
///
/// * `spot_price` - The current price of the underlying asset.
/// * `strike_price` - The strike price of the option.
/// * `risk_free_rate` - The risk-free interest rate.
/// * `volatility` - The volatility of the underlying asset.
/// * `time_to_maturity` - The time to maturity of the option (in years).
/// * `num_simulations` - The number of simulations to run for the Monte Carlo method.
/// * `num_steps` - The number of steps in each simulation for the Monte Carlo method.
///
/// # Returns
///
/// The rho of the American call option.
///
/// # Example
///
/// ```rust
/// use monte_carlo_options_rs::{calculate_rho, monte_carlo_american_call};
///
/// let spot_price = 100.0;
/// let strike_price = 110.0;
/// let risk_free_rate = 0.05;
/// let volatility = 0.2;
/// let time_to_maturity = 1.0;
/// let num_simulations = 100_000;
/// let num_steps = 252;
///
/// let rho = calculate_rho(
///     spot_price,
///     strike_price,
///     risk_free_rate,
///     volatility,
///     time_to_maturity,
///     num_simulations,
///     num_steps,
/// ).unwrap();
///
/// println!("Rho of the American call option: {:.4}", rho);
/// ```
#[pyfunction]
fn calculate_rho(
    spot_price: f64,
    strike_price: f64,
    risk_free_rate: f64,
    volatility: f64,
    time_to_maturity: f64,
    num_simulations: usize,
    num_steps: usize,
) -> PyResult<f64> {
    let h = 0.01;
    let price_up = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate * (1.0 + h),
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    let price_down = monte_carlo_american_call(
        spot_price,
        strike_price,
        risk_free_rate * (1.0 - h),
        volatility,
        time_to_maturity,
        num_simulations,
        num_steps,
    )?;
    Ok((price_up - price_down) / (2.0 * h))
}

/// Initializes the `libmonte_carlo_options_rs` Python module.
///
/// This function is called when the Python module is imported and adds the Rust functions
/// to the module so they can be accessed from Python.
///
/// # Arguments
///
/// * `_py` - The Python interpreter.
/// * `m` - The Python module to add the functions to.
///
/// # Returns
///
/// `Ok(())` if the module was initialized successfully, or an error if there was a problem.
///
/// # Example
///
/// ```python
/// import libmonte_carlo_options_rs
///
/// # Use the functions from the module
/// price = libmonte_carlo_options_rs.monte_carlo_american_call(...)
/// delta = libmonte_carlo_options_rs.calculate_delta(...)
/// gamma = libmonte_carlo_options_rs.calculate_gamma(...)
/// vega = libmonte_carlo_options_rs.calculate_vega(...)
/// theta = libmonte_carlo_options_rs.calculate_theta(...)
/// rho = libmonte_carlo_options_rs.calculate_rho(...)
/// ```
#[pymodule]
fn libmonte_carlo_options_rs(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(monte_carlo_american_call, m)?)?;
    m.add_function(wrap_pyfunction!(monte_carlo_american_put, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_delta, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_gamma, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_vega, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_theta, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_rho, m)?)?;
    Ok(())
}