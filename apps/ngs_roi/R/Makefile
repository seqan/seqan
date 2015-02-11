VERSION=0.2
PKG=ngsroi_${VERSION}.tar.gz

all: ${PKG} check

package: ${PKG}

${PKG}:
	R CMD build ngsroi

check: ${PKG}
	R CMD check ngsroi
