# Optimization Library

This C library provides implementations of various optimization algorithms for finding roots and minimizing functions.

## Features

- **Brent's Method for Root Finding**: Efficiently finds roots of one-dimensional functions.
- **Brent's Method for Minimization**: Minimizes one-dimensional functions.
- **Nelder-Mead Method**: Minimizes multidimensional functions without requiring derivatives.

## Installation

1. Clone the repository:
   ```
   git clone https://github.com/agentlans/optimization-library.git
   ```
2. Navigate to the project directory:
   ```
   cd optimization-library
   ```
3. Compile the library using the provided Makefile:
   ```
   make
   ```

## Usage

Include the `opt.h` header in your C program:

```c
#include "opt.h"
```

### Brent's Method for Root Finding

```c
double root = OPT_BrentRoot(function, lower_bound, upper_bound, params, tolerance, max_iterations, &error);
```

### Brent's Method for Minimization

```c
double min_value = OPT_BrentMinimize(function, lower_bound, upper_bound, params, tolerance, max_iterations, &x_min, &error);
```

### Nelder-Mead Method

```c
OPT_NelderMead(function, x_min, dimensions, params, tolerance, max_iterations, initial_step, &error);
```

## Error Handling

The library uses an `OPT_Error` enum to indicate the status of operations. Always check the `error` parameter after calling a function to ensure successful execution.

## Example

See `example.c` for sample usage of the library functions.

## Building and Running the Example

```
make example
./example
```

## Contact

For any questions or concerns, please open an issue on the GitHub repository.

## License

Copyright :copyright: 2024 by Alan Tseng

[MIT License](LICENSE).
