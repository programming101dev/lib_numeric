/*
 * Copyright 2026 D'Arcy Smith.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "p101_math/p101_math.h"
#include <fenv.h>
#include <math.h>
#include <p101_c/p101_fenv.h>
#include <p101_c/p101_math.h>
#include <p101_env/wrapper.h>

static int bessel_error_code(const struct p101_env *env);
static int bessel_y_error_code(const struct p101_env *env, double x);
static int prepare_bessel_call(const struct p101_env *env, struct p101_error *err);

static int bessel_error_code(const struct p101_env *env)
{
    int p101_single_result_;
    int exceptions;

    if(errno != 0)
    {
        p101_single_result_ = errno;
        goto p101_single_exit_;
    }

    exceptions = p101_fetestexcept(env, FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW);

    if((exceptions & FE_INVALID) != 0)
    {
        p101_single_result_ = EDOM;
        goto p101_single_exit_;
    }

    if((exceptions & (FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW)) != 0)
    {
        p101_single_result_ = ERANGE;
        goto p101_single_exit_;
    }

    p101_single_result_ = 0;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

static int bessel_y_error_code(const struct p101_env *env, double x)
{
    int p101_single_result_;
    int error_code;

    error_code = bessel_error_code(env);
    if(error_code != 0)
    {
        p101_single_result_ = error_code;
        goto p101_single_exit_;
    }

    /*
     * Some supported libm implementations return the required exceptional
     * value without setting errno or a floating-point exception. The XSI
     * domain itself is still unambiguous.
     */
    if(x < 0.0)
    {
        p101_single_result_ = EDOM;
        goto p101_single_exit_;
    }

    if(fpclassify(x) == FP_ZERO)
    {
        p101_single_result_ = ERANGE;
        goto p101_single_exit_;
    }

    p101_single_result_ = 0;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

static int prepare_bessel_call(const struct p101_env *env, struct p101_error *err)
{
    int result;

    result = p101_feclearexcept(env, err, FE_ALL_EXCEPT);
    if(result == 0)
    {
        /*
         * p101_feclearexcept() is itself instrumented.  Keep any errno used
         * by that instrumentation from becoming the native math operation's
         * apparent failure.
         */
        errno = 0;
    }
    return result;
}

double p101_j0(const struct p101_env *env, struct p101_error *err, double x)
{
    double p101_single_result_;
    int    prepare_result;
    int    error_code;
    double ret_val;

    P101_TRACE(env);
    P101_WRAPPER_FAULT_RETURN(env, err, ret_val, p101_nan(env, ""));
    prepare_result = prepare_bessel_call(env, err);
    if(prepare_result != 0)
    {
        P101_WRAPPER_DONE(env);
        p101_single_result_ = p101_nan(env, "");
        goto p101_single_exit_;
    }
    ret_val    = j0(x);
    error_code = bessel_error_code(env);

    if(error_code != 0)
    {
        P101_ERROR_RAISE_ERRNO(err, error_code);
    }

    P101_TRACE_EXIT(env);
    p101_single_result_ = ret_val;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

double p101_j1(const struct p101_env *env, struct p101_error *err, double x)
{
    double p101_single_result_;
    int    prepare_result;
    int    error_code;
    double ret_val;

    P101_TRACE(env);
    P101_WRAPPER_FAULT_RETURN(env, err, ret_val, p101_nan(env, ""));
    prepare_result = prepare_bessel_call(env, err);
    if(prepare_result != 0)
    {
        P101_WRAPPER_DONE(env);
        p101_single_result_ = p101_nan(env, "");
        goto p101_single_exit_;
    }
    ret_val    = j1(x);
    error_code = bessel_error_code(env);

    if(error_code != 0)
    {
        P101_ERROR_RAISE_ERRNO(err, error_code);
    }

    P101_TRACE_EXIT(env);
    p101_single_result_ = ret_val;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

double p101_jn(const struct p101_env *env, struct p101_error *err, int n, double x)
{
    double p101_single_result_;
    int    prepare_result;
    int    error_code;
    double ret_val;

    P101_TRACE(env);
    P101_WRAPPER_FAULT_RETURN(env, err, ret_val, p101_nan(env, ""));
    prepare_result = prepare_bessel_call(env, err);
    if(prepare_result != 0)
    {
        P101_WRAPPER_DONE(env);
        p101_single_result_ = p101_nan(env, "");
        goto p101_single_exit_;
    }
    ret_val    = jn(n, x);
    error_code = bessel_error_code(env);

    if(error_code != 0)
    {
        P101_ERROR_RAISE_ERRNO(err, error_code);
    }

    P101_TRACE_EXIT(env);
    p101_single_result_ = ret_val;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

double p101_y0(const struct p101_env *env, struct p101_error *err, double x)
{
    double p101_single_result_;
    int    prepare_result;
    int    error_code;
    double ret_val;

    P101_TRACE(env);
    P101_WRAPPER_FAULT_RETURN(env, err, ret_val, p101_nan(env, ""));
    prepare_result = prepare_bessel_call(env, err);
    if(prepare_result != 0)
    {
        P101_WRAPPER_DONE(env);
        p101_single_result_ = p101_nan(env, "");
        goto p101_single_exit_;
    }
    ret_val    = y0(x);
    error_code = bessel_y_error_code(env, x);

    if(error_code != 0)
    {
        P101_ERROR_RAISE_ERRNO(err, error_code);
    }

    P101_TRACE_EXIT(env);
    p101_single_result_ = ret_val;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

double p101_y1(const struct p101_env *env, struct p101_error *err, double x)
{
    double p101_single_result_;
    int    prepare_result;
    int    error_code;
    double ret_val;

    P101_TRACE(env);
    P101_WRAPPER_FAULT_RETURN(env, err, ret_val, p101_nan(env, ""));
    prepare_result = prepare_bessel_call(env, err);
    if(prepare_result != 0)
    {
        P101_WRAPPER_DONE(env);
        p101_single_result_ = p101_nan(env, "");
        goto p101_single_exit_;
    }
    ret_val    = y1(x);
    error_code = bessel_y_error_code(env, x);

    if(error_code != 0)
    {
        P101_ERROR_RAISE_ERRNO(err, error_code);
    }

    P101_TRACE_EXIT(env);
    p101_single_result_ = ret_val;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}

double p101_yn(const struct p101_env *env, struct p101_error *err, int n, double x)
{
    double p101_single_result_;
    int    prepare_result;
    int    error_code;
    double ret_val;

    P101_TRACE(env);
    P101_WRAPPER_FAULT_RETURN(env, err, ret_val, p101_nan(env, ""));
    prepare_result = prepare_bessel_call(env, err);
    if(prepare_result != 0)
    {
        P101_WRAPPER_DONE(env);
        p101_single_result_ = p101_nan(env, "");
        goto p101_single_exit_;
    }
    ret_val    = yn(n, x);
    error_code = bessel_y_error_code(env, x);

    if(error_code != 0)
    {
        P101_ERROR_RAISE_ERRNO(err, error_code);
    }

    P101_TRACE_EXIT(env);
    p101_single_result_ = ret_val;
    goto p101_single_exit_;

p101_single_exit_:
    return p101_single_result_;
}
