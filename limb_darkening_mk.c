#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double calc_mu(gsl_rng *);
double calc_initial_mu(gsl_rng *);
double calc_tau(gsl_rng *);
int do_random_walk(int n_photons, double tau_max, int n_mu, int mu_array[]);


double calc_mu(gsl_rng *rng) {
  return 2.0 * ( gsl_rng_uniform(rng) - 0.5 );
}

double calc_initial_mu(gsl_rng *rng) {
  return sqrt(gsl_rng_uniform(rng));
}

double calc_tau(gsl_rng *rng) {
  return gsl_ran_exponential(rng, 1.0);
}


int do_random_walk(int n_photons, double tau_max, int n_mu, int
mu_array[]) {
  int i, j;
  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default); // See the GSL manual
for ways of controlling the RNG through environment variables

  for (i = 0; i < n_mu; ++i) mu_array[i] = 0;

  for (i = 0; i < n_photons; ++i) {
    double mu; // cos theta
    double tau; // optical depth of the photon

    // start at the deepest layer
    tau = tau_max;
    // compute the initial emission angle
    mu = calc_initial_mu(rng);
    // first step
    tau = tau - calc_tau(rng) * mu;

    while (tau > 0.0) {
      if (tau > tau_max) {
    // If we left the photosphere through the inner boundary, emit another
photon.
    mu = calc_initial_mu(rng);
    tau = tau - calc_tau(rng) * mu;
      }
      else {
    mu = calc_mu(rng);
    tau = tau - calc_tau(rng) * mu;
      }
    }

    // find the correct bin; the final mu is in the interval (0,1)
    j = floor(mu * n_mu);
    mu_array[j]++;

  }

  gsl_rng_free(rng);

  return 0;
}

int main() {
  int n_photons; //  number of photons
  double tau_max; // optical detph at which the simulation starts
  const int n_mu = 100; // number of angular bins
  int *mu_array; // angular bins
  int i;
  char* filename = "output";
  FILE *file;

  gsl_rng_env_setup();

  mu_array = malloc(n_mu*sizeof(int));
  if (!mu_array) {
    perror("Could not allocate storage for the results.");
    exit(EXIT_FAILURE);
  }

  n_photons = 1e5;
  tau_max = 10.0;
  do_random_walk(n_photons, tau_max, n_mu, mu_array);

  file = fopen(filename, "w");
  if (!file) {
    fprintf(stderr, "Error opening `%s' for writing. Writing to stdout
instead.", filename);
    file = stdout;
  }

  for (i = 0; i < n_mu; ++i) fprintf(file, "%g %g\n", (i +
0.5)/((double) n_mu), (double) mu_array[i] / (double) mu_array[n_mu-1]);

  fclose(file);

  free(mu_array);
  return 0;
}

