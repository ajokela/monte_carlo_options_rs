# Monte Carlo Options Pricing in Rust

This repository contains a Rust library for pricing American call and put options using the Monte Carlo simulation method. The library also provides functions to calculate various option Greeks, such as delta, gamma, vega, theta, and rho, using the finite difference method.

## Features

- Price American call options using Monte Carlo simulation
- Price American put options using Monte Carlo simulation
- Calculate option Greeks: delta, gamma, vega, theta, and rho
- Expose the library functions to Python using PyO3

## Installation

To use this library in your Rust project, add the following to your `Cargo.toml` file:

```toml
[dependencies]
monte_carlo_options_rs = { git = "https://github.com/ajokela/monte_carlo_options_rs.git" }
```

## Usage

### Rust

```rust
use monte_carlo_options_rs::{monte_carlo_american_call, calculate_delta};

let spot_price = 100.0;
let strike_price = 110.0;
let risk_free_rate = 0.05;
let volatility = 0.2;
let time_to_maturity = 1.0;
let num_simulations = 100_000;
let num_steps = 252;

let price = monte_carlo_american_call(
    spot_price,
    strike_price,
    risk_free_rate,
    volatility,
    time_to_maturity,
    num_simulations,
    num_steps,
).unwrap();

let delta = calculate_delta(
    spot_price,
    strike_price,
    risk_free_rate,
    volatility,
    time_to_maturity,
    num_simulations,
    num_steps,
).unwrap();

println!("American call option price: {:.2}", price);
println!("Delta of the American call option: {:.4}", delta);
```

### Python

```python
import libmonte_carlo_options_rs

spot_price = 100.0
strike_price = 110.0
risk_free_rate = 0.05
volatility = 0.2
time_to_maturity = 1.0
num_simulations = 100_000
num_steps = 252

price = libmonte_carlo_options_rs.monte_carlo_american_call(
    spot_price,
    strike_price,
    risk_free_rate,
    volatility,
    time_to_maturity,
    num_simulations,
    num_steps,
)

delta = libmonte_carlo_options_rs.calculate_delta(
    spot_price,
    strike_price,
    risk_free_rate,
    volatility,
    time_to_maturity,
    num_simulations,
    num_steps,
)

print(f"American call option price: {price:.2f}")
print(f"Delta of the American call option: {delta:.4f}")
```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

This project is licensed under the terms of the 3-Clause BSD license. See [LICENSE](LICENSE) for more details.