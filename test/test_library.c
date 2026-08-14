#include <p101_env/env.h>
#include <p101_error/error.h>
#include <p101_math/p101_math.h>
#include <p101_random/p101_stdlib.h>
#include <stdlib.h>

int main(void)
{
    struct p101_error *err;
    struct p101_env   *env;
    int                exit_status;

    exit_status = EXIT_FAILURE;
    err         = p101_error_create(false);
    if(err != NULL)
    {
        env = p101_env_create(err, NULL);
        if(env != NULL)
        {
            exit_status = EXIT_SUCCESS;
            p101_env_destroy(env);
        }
        p101_error_destroy(err);
    }

    return exit_status;
}
