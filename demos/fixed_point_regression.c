#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "STFT.h"
#include "fixed_math.h"
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef BTT_REGRESSION_TOLERANCE
#define BTT_REGRESSION_TOLERANCE 0.05f
#endif

#ifndef BTT_REGRESSION_WINDOW
#define BTT_REGRESSION_WINDOW 512
#endif

#ifndef BTT_REGRESSION_OVERLAP
#define BTT_REGRESSION_OVERLAP 4
#endif

#ifndef BTT_REGRESSION_FRAMES_MAX
#define BTT_REGRESSION_FRAMES_MAX 256
#endif

static void generate_sinusoid(float* dst, int len, double freq_hz, double sample_rate)
{
  double phase = 0.0;
  double step = (2.0 * M_PI * freq_hz) / sample_rate;
  for(int i = 0; i < len; ++i)
    {
      dst[i] = (float)sin(phase);
      phase += step;
    }
}

static void generate_impulse(float* dst, int len)
{
  memset(dst, 0, (size_t)len * sizeof(*dst));
  if(len > 0)
    {
      dst[0] = 1.0f;
    }
}

static void generate_percussive_burst(float* dst, int len, int burst_len)
{
  memset(dst, 0, (size_t)len * sizeof(*dst));
  if(burst_len > len)
    {
      burst_len = len;
    }
  for(int i = 0; i < burst_len; ++i)
    {
      float env = expf(-0.01f * (float)i);
      float noise = (float)rand() / (float)RAND_MAX * 2.0f - 1.0f;
      dst[i] = env * noise;
    }
}

struct baseline_frames
{
  dft_sample_t* frames;
  int           frame_count;
  int           fft_N;
};

struct compare_stats
{
  const struct baseline_frames* baseline;
  int           frame_index;
  double        mse_accum;
  double        max_error;
  int           samples;
  float         tolerance;
  int           failures;
};

static void capture_baseline(void* ctx, dft_sample_t* real, int N)
{
  struct baseline_frames* b = (struct baseline_frames*)ctx;
  if(b->frame_count >= BTT_REGRESSION_FRAMES_MAX)
    {
      return;
    }
  memcpy(&b->frames[b->frame_count * N], real, (size_t)N * sizeof(*real));
  b->frame_count++;
  b->fft_N = N;
}

static void compare_against_baseline(void* ctx, dft_sample_t* real, int N)
{
  struct compare_stats* stats = (struct compare_stats*)ctx;
  if(stats->frame_index >= stats->baseline->frame_count || stats->frame_index >= BTT_REGRESSION_FRAMES_MAX)
    {
      return;
    }
  const dft_sample_t* ref = &stats->baseline->frames[stats->frame_index * N];
  double ref_energy = 0.0;
  double cand_energy = 0.0;
  double dot = 0.0;
  for(int i = 0; i < N; ++i)
    {
      ref_energy  += (double)ref[i]  * (double)ref[i];
      cand_energy += (double)real[i] * (double)real[i];
      dot         += (double)ref[i]  * (double)real[i];
    }
  double gain = 1.0;
  if(cand_energy > 0.0)
    {
      gain = sqrt(ref_energy / cand_energy);
    }
  double frame_mse = 0.0;
  double aligned_energy = 0.0;
  double dot_aligned = 0.0;
  for(int i = 0; i < N; ++i)
    {
      double aligned = (double)real[i] * gain;
      double diff = aligned - (double)ref[i];
      frame_mse += diff * diff;
      aligned_energy += aligned * aligned;
      dot_aligned += aligned * (double)ref[i];
      if(fabs(diff) > stats->max_error)
        {
          stats->max_error = fabs(diff);
        }
    }
  stats->mse_accum += frame_mse;
  double corr = 1.0;
  if(ref_energy > 0.0 && aligned_energy > 0.0)
    {
      corr = dot_aligned / sqrt(ref_energy * aligned_energy);
    }
  double frame_rmse = sqrt(frame_mse / (double)N);
  if(frame_rmse > stats->tolerance && corr < (1.0 - stats->tolerance))
    {
      stats->failures++;
    }
  stats->samples += N;
  stats->frame_index++;
}

static void noop_onprocess(void* ctx, dft_sample_t* real, int N)
{
  (void)ctx;
  (void)real;
  (void)N;
}

static int run_regression(const char* name, float* input, int len, double sample_rate, float tolerance, double* out_mse, double* out_max_err)
{
  int hop = BTT_REGRESSION_WINDOW / BTT_REGRESSION_OVERLAP;
  int max_frames = (len + hop - 1) / hop;
  if(max_frames > BTT_REGRESSION_FRAMES_MAX)
    {
      max_frames = BTT_REGRESSION_FRAMES_MAX;
    }

  struct baseline_frames baseline;
  baseline.frames = (dft_sample_t*)calloc((size_t)max_frames * BTT_REGRESSION_WINDOW, sizeof(*baseline.frames));
  baseline.frame_count = 0;
  baseline.fft_N = 0;

  struct compare_stats stats;
  stats.baseline  = &baseline;
  stats.frame_index = 0;
  stats.mse_accum = 0.0;
  stats.max_error = 0.0;
  stats.samples   = 0;
  stats.tolerance = tolerance;
  stats.failures  = 0;

  STFT* stft_float = stft_new(BTT_REGRESSION_WINDOW, BTT_REGRESSION_OVERLAP, 0);
  STFT* stft_fixed = stft_new(BTT_REGRESSION_WINDOW, BTT_REGRESSION_OVERLAP, 0);
  if(stft_float == NULL || stft_fixed == NULL || baseline.frames == NULL)
    {
      fprintf(stderr, "%s: allocation failed\n", name);
      free(baseline.frames);
      stft_destroy(stft_float);
      stft_destroy(stft_fixed);
      return -1;
    }

  stft_set_use_fixed_point(stft_fixed, 1);

  float* float_input = (float*)malloc((size_t)len * sizeof(*float_input));
  float* fixed_input = (float*)malloc((size_t)len * sizeof(*fixed_input));
  if(float_input == NULL || fixed_input == NULL)
    {
      fprintf(stderr, "%s: allocation failed\n", name);
      free(baseline.frames);
      stft_destroy(stft_float);
      stft_destroy(stft_fixed);
      free(float_input);
      free(fixed_input);
      return -1;
    }

  memcpy(float_input, input, (size_t)len * sizeof(*float_input));
  memcpy(fixed_input, input, (size_t)len * sizeof(*fixed_input));

  stft_process(stft_float, float_input, len, capture_baseline, &baseline);
  stft_process(stft_fixed, fixed_input, len, compare_against_baseline, &stats);

  *out_mse     = (stats.samples > 0) ? (stats.mse_accum / (double)stats.samples) : 0.0;
  *out_max_err = stats.max_error;

  printf("%s: sr=%.1fHz frames=%d mse=%g max_err=%g tolerance=%g failures=%d\n",
         name, sample_rate, stats.frame_index, *out_mse, *out_max_err, tolerance, stats.failures);

  free(baseline.frames);
  stft_destroy(stft_float);
  stft_destroy(stft_fixed);
  free(float_input);
  free(fixed_input);

  return stats.failures == 0 ? 0 : 1;
}

static void profile_stft(float* input, int len, double sample_rate)
{
  STFT* stft = stft_new(BTT_REGRESSION_WINDOW, BTT_REGRESSION_OVERLAP, 0);
  if(stft == NULL)
    {
      return;
    }
  stft_set_use_fixed_point(stft, 1);
  clock_t start = clock();
  stft_process(stft, input, len, noop_onprocess, NULL);
  clock_t end = clock();
  double seconds = (double)(end - start) / (double)CLOCKS_PER_SEC;
  double samples_per_second = (seconds > 0.0) ? (double)len / seconds : 0.0;
  double cycles_per_sample = seconds > 0.0 ? ((double)CLOCKS_PER_SEC / samples_per_second) : 0.0;
  printf("profile: len=%d sr=%.2fHz time=%.6fs samples/s=%.2f host-cycles/sample=%.2f\n",
         len, sample_rate, seconds, samples_per_second, cycles_per_sample);
  stft_destroy(stft);
}

int main(void)
{
  const int len = 4096;
  float* buffer = (float*)malloc((size_t)len * sizeof(*buffer));
  if(buffer == NULL)
    {
      return 1;
    }

  double mse, max_err;
  int rc = 0;

  generate_sinusoid(buffer, len, 440.0, 44100.0);
  rc |= run_regression("sinusoid_44k1", buffer, len, 44100.0, BTT_REGRESSION_TOLERANCE, &mse, &max_err);

  generate_sinusoid(buffer, len, 220.5, 22050.0);
  rc |= run_regression("sinusoid_22k05", buffer, len, 22050.0, BTT_REGRESSION_TOLERANCE, &mse, &max_err);

  generate_impulse(buffer, len);
  rc |= run_regression("impulse", buffer, len, 44100.0, BTT_REGRESSION_TOLERANCE, &mse, &max_err);

  generate_percussive_burst(buffer, len, 64);
  rc |= run_regression("percussive_burst", buffer, len, 44100.0, BTT_REGRESSION_TOLERANCE, &mse, &max_err);

  profile_stft(buffer, len, 22050.0);
  profile_stft(buffer, len, 44100.0);

  free(buffer);

  if(rc != 0)
    {
      fprintf(stderr, "Fixed-point regression: FAIL (see failures above)\n");
    }
  else
    {
      printf("Fixed-point regression: PASS\n");
    }

  return rc;
}

