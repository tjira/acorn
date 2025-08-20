#!/bin/bash

CORES=$(nproc --all); DIRNAME=$(dirname "$(realpath $0)"); TARGETS=("x86_64-linux")

CMAKE_GENERAL=(
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_C_COMPILER="$DIRNAME/../zigcc"
    -DCMAKE_CXX_COMPILER="$DIRNAME/../zigcpp"
)

CMAKE_LIBINT=(
    -DLIBINT_USE_BUNDLED_BOOST=True
)

CMAKE_OMP=(
    -DLIBOMP_ENABLE_SHARED=False
    -DLIBOMP_OMPD_SUPPORT=False
    -DLIBOMP_OMPT_SUPPORT=False
)

CMAKE_RAYLIB=(
    -DPLATFORM=Desktop
)

CONFIGURE_GENERAL=(
    CC="$DIRNAME/../zigcc"
    --enable-static
    --disable-shared
)

MAKE_OPENBLAS=(
    AR="$DIRNAME/../zigar"
    CC="$DIRNAME/../zigcc"
    DYNAMIC_ARCH=1
    HOSTCC=gcc
    NO_SHARED=1
    NOFORTRAN=1
    NUM_THREADS=128
    RANLIB="$DIRNAME/../zigranlib"
    USE_OPENMP=1
)

cmake_install() {
    cmake --build . --parallel $CORES --verbose && cmake --install .
}

tar -xf external-x86_64-linux/zig.tar.xz -C zig-bin --strip-components=1
tar -xf external-x86_64-linux/zls.tar.xz -C zig-bin --strip-components=0

for TARGET in "${TARGETS[@]}"; do

    export    ACLOCAL_PATH="$DIRNAME/../external-$TARGET/share/aclocal:$ACLOCAL_PATH"
    export  C_INCLUDE_PATH="$DIRNAME/../external-$TARGET/include:$C_INCLUDE_PATH"
    export PKG_CONFIG_PATH="$DIRNAME/../external-$TARGET/share/pkgconfig:$PKG_CONFIG_PATH"
    export PKG_CONFIG_PATH="$DIRNAME/../external-$TARGET/lib/pkgconfig:$PKG_CONFIG_PATH"

    CMAKE_GENERAL+=(
        -DCMAKE_INSTALL_PREFIX="$DIRNAME/../external-$TARGET"
        -DCMAKE_PREFIX_PATH="$DIRNAME/../external-$TARGET"
    )

    CONFIGURE_GENERAL+=(
        --prefix="$DIRNAME/../external-$TARGET"
    )

    echo "$DIRNAME/../zig-bin/zig ar                          \"\$@\"" > zigar     && chmod +x     zigar
    echo "$DIRNAME/../zig-bin/zig cc     --target=$TARGET-gnu \"\$@\"" > zigcc     && chmod +x     zigcc
    echo "$DIRNAME/../zig-bin/zig c++    --target=$TARGET-gnu \"\$@\"" > zigcpp    && chmod +x    zigcpp
    echo "$DIRNAME/../zig-bin/zig ranlib                      \"\$@\"" > zigranlib && chmod +x zigranlib

    tar -xf  external-$TARGET/llvm.tar.xz       && mv llvm*            llvm
    tar -xzf external-$TARGET/eigen.tar.gz      && mv eigen*          eigen
    tar -xzf external-$TARGET/libint.tar.gz     && mv libint*        libint
    tar -xzf external-$TARGET/mesa.tar.gz       && mv mesa*            mesa
    tar -xzf external-$TARGET/openblas.tar.gz   && mv OpenBLAS*    openblas
    tar -xzf external-$TARGET/raylib.tar.gz     && mv raylib*        raylib
    tar -xzf external-$TARGET/x11.tar.gz        && mv libx11*           x11
    tar -xzf external-$TARGET/xau.tar.gz        && mv libxau*           xau
    tar -xzf external-$TARGET/xcb.tar.gz        && mv libxcb*           xcb
    tar -xzf external-$TARGET/xcbproto.tar.gz   && mv xcbproto*    xcbproto
    tar -xzf external-$TARGET/xcursor.tar.gz    && mv libxcursor*   xcursor
    tar -xzf external-$TARGET/xext.tar.gz       && mv libxext*         xext
    tar -xzf external-$TARGET/xfixes.tar.gz     && mv libxfixes*     xfixes
    tar -xzf external-$TARGET/xi.tar.gz         && mv libxi*             xi
    tar -xzf external-$TARGET/xinerama.tar.gz   && mv libxinerama* xinerama
    tar -xzf external-$TARGET/xinput.tar.gz     && mv xinput*        xinput
    tar -xzf external-$TARGET/xorgmacros.tar.gz && mv macros*    xorgmacros
    tar -xzf external-$TARGET/xorgproto.tar.gz  && mv xorgproto*  xorgproto
    tar -xzf external-$TARGET/xrandr.tar.gz     && mv libxrandr*     xrandr
    tar -xzf external-$TARGET/xrender.tar.gz    && mv libxrender*   xrender
    tar -xzf external-$TARGET/xtrans.tar.gz     && mv libxtrans*     xtrans

    cd xorgmacros && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xorgproto  && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xtrans     && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xau        && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xcbproto   && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xcb        && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd x11        && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xext       && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xrender    && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xrandr     && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xinerama   && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xfixes     && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xcursor    && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xi         && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xinput     && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..

    # cd mesa && mkdir build && meson setup build -Dprefix="$DIRNAME/../external-$TARGET"

    # cd llvm   && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_OMP[@]}"    ../openmp && cmake_install && cd ../..
    # cd eigen  && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}"                      ..        && cmake_install && cd ../..
    # cd libint && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_LIBINT[@]}" ..        && cmake_install && cd ../..
    cd raylib && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_RAYLIB[@]}" ..        && cmake_install && cd ../..

    # cd openblas && make "${MAKE_OPENBLAS[@]}" -j $CORES libs shared && make NO_SHARED=1 PREFIX="$DIRNAME/../external-$TARGET" install && cd ..

    # rm -rf eigen libint llvm mesa openblas raylib x11 xau xcb xcbproto xext xfixes xi xinerama xinput xorgmacros xorgproto xrender xrandr xtrans zigar zigcc zigcpp zigranlib

done
