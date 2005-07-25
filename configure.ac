AC_PREREQ(2.50)
AC_INIT([DESCRIPTION])

: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CXX=`${R_HOME}/bin/R CMD config CXX`

AC_PROG_CXX

if test "${GXX}" = yes; then
  gxx_version=`${CXX} --version`
  case ${gxx_version} in
    2.96*)
      AC_MSG_WARN([g++ 2.96 cannot reliably be used with this package.])
      AC_MSG_ERROR([Please use a newer version of g++ or a different C++ compiler.])   
   ;;
  esac
fi
