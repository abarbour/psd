# multitaper

Version: 1.0-14

## In both

*   checking whether package ‘multitaper’ can be installed ... ERROR
    ```
    Installation failed.
    See ‘/Users/abl/survey.processing/development/R/packages/psd/revdep/checks.noindex/multitaper/new/multitaper.Rcheck/00install.out’ for details.
    ```

## Installation

### Devel

```
* installing *source* package ‘multitaper’ ...
** package ‘multitaper’ successfully unpacked and MD5 sums checked
** libs
gfortran   -fPIC  -g -O2  -c djt.f -o djt.o
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
gfortran  -fPIC -Wall -g -O2  -c  dpss.f90 -o dpss.o
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c multitaper_init.c -o multitaper_init.o
gfortran   -fPIC  -g -O2  -c sine.f -o sine.o
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
gfortran -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o multitaper.so djt.o dpss.o multitaper_init.o sine.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0'
ld: warning: directory not found for option '-L/usr/local/gfortran/lib'
ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0'
ld: warning: directory not found for option '-L/usr/local/gfortran/lib'
ld: library not found for -ldylib1.o
collect2: error: ld returned 1 exit status
make: *** [multitaper.so] Error 1
ERROR: compilation failed for package ‘multitaper’
* removing ‘/Users/abl/survey.processing/development/R/packages/psd/revdep/checks.noindex/multitaper/new/multitaper.Rcheck/multitaper’

```
### CRAN

```
* installing *source* package ‘multitaper’ ...
** package ‘multitaper’ successfully unpacked and MD5 sums checked
** libs
gfortran   -fPIC  -g -O2  -c djt.f -o djt.o
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
gfortran  -fPIC -Wall -g -O2  -c  dpss.f90 -o dpss.o
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
clang -I"/Library/Frameworks/R.framework/Resources/include" -DNDEBUG   -I/usr/local/include   -fPIC  -Wall -g -O2  -c multitaper_init.c -o multitaper_init.o
gfortran   -fPIC  -g -O2  -c sine.f -o sine.o
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
gfortran -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib -o multitaper.so djt.o dpss.o multitaper_init.o sine.o -L/Library/Frameworks/R.framework/Resources/lib -lRlapack -L/Library/Frameworks/R.framework/Resources/lib -lRblas -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm -L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0 -L/usr/local/gfortran/lib -lgfortran -lquadmath -lm -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
gfortran: warning: couldn’t understand kern.osversion ‘18.2.0
ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0'
ld: warning: directory not found for option '-L/usr/local/gfortran/lib'
ld: warning: directory not found for option '-L/usr/local/gfortran/lib/gcc/x86_64-apple-darwin15/6.1.0'
ld: warning: directory not found for option '-L/usr/local/gfortran/lib'
ld: library not found for -ldylib1.o
collect2: error: ld returned 1 exit status
make: *** [multitaper.so] Error 1
ERROR: compilation failed for package ‘multitaper’
* removing ‘/Users/abl/survey.processing/development/R/packages/psd/revdep/checks.noindex/multitaper/old/multitaper.Rcheck/multitaper’

```
