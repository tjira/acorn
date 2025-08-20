#!/bin/bash

CORES=$(nproc --all); DIRNAME=$(dirname "$(realpath $0)"); TARGETS=("x86_64-linux")

CMAKE_GENERAL=(
    -DBUILD_SHARED_LIBS=OFF
    -DCMAKE_BUILD_TYPE=Release
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
    --enable-static
    --disable-shared
)

MAKE_OPENBLAS=(
    DYNAMIC_ARCH=1
    HOSTCC=gcc
    NO_SHARED=1
    NOFORTRAN=1
    NUM_THREADS=128
    USE_OPENMP=1
)

MESON_GENERAL=(
    -Dbuildtype=release
    -Ddefault_library=static
)

MESON_DRM=(
    -Dintel=enabled
)

MESON_MESA=(
    -Dosmesa=true
    -Damber=true
    -Dglvnd=false
    -Dgallium-drivers=swrast
    -Dplatforms=x11
    -Dglx=disabled
    -Dshared-glapi=disabled
    -Ddri-drivers=
    -Degl=disabled
    -Dgbm=disabled
    -Dvulkan-drivers=
)

cmake_install() {
    cmake --build . --parallel $CORES --verbose && cmake --install .
}

tar -xf external-x86_64-linux/zig.tar.xz -C zig-bin --strip-components=1
tar -xf external-x86_64-linux/zls.tar.xz -C zig-bin --strip-components=0

export CC="$DIRNAME/../zigcc"; export CXX="$DIRNAME/../zigcpp"; export AR="$DIRNAME/../zigar"; export RANLIB="$DIRNAME/../zigranlib"

for TARGET in "${TARGETS[@]}"; do

    export    ACLOCAL_PATH="$DIRNAME/../external-$TARGET/share/aclocal:$ACLOCAL_PATH"
    export  C_INCLUDE_PATH="$DIRNAME/../external-$TARGET/include:$C_INCLUDE_PATH"

    export PKG_CONFIG_PATH="$DIRNAME/../external-$TARGET/share/pkgconfig:$PKG_CONFIG_PATH"
    export PKG_CONFIG_PATH="$DIRNAME/../external-$TARGET/lib/pkgconfig:$PKG_CONFIG_PATH"
    export PKG_CONFIG_PATH="$DIRNAME/../external-$TARGET/lib/$TARGET-gnu/pkgconfig:$PKG_CONFIG_PATH"

    CMAKE_GENERAL+=(
        -DCMAKE_INSTALL_PREFIX="$DIRNAME/../external-$TARGET"
        -DCMAKE_PREFIX_PATH="$DIRNAME/../external-$TARGET"
    )

    CONFIGURE_GENERAL+=(
        --prefix="$DIRNAME/../external-$TARGET"
    )

    MESON_GENERAL+=(
        -Dprefix="$DIRNAME/../external-$TARGET"
    )

    echo -e "#!/usr/bin/env bash\n\n$DIRNAME/../zig-bin/zig ar                          \"\$@\"" > zigar     && chmod +x     zigar
    echo -e "#!/usr/bin/env bash\n\n$DIRNAME/../zig-bin/zig cc     --target=$TARGET-gnu \"\$@\"" > zigcc     && chmod +x     zigcc
    echo -e "#!/usr/bin/env bash\n\n$DIRNAME/../zig-bin/zig c++    --target=$TARGET-gnu \"\$@\"" > zigcpp    && chmod +x    zigcpp
    echo -e "#!/usr/bin/env bash\n\n$DIRNAME/../zig-bin/zig ranlib                      \"\$@\"" > zigranlib && chmod +x zigranlib

    tar -xf  external-$TARGET/llvm.tar.xz       && mv llvm*            llvm
    tar -xzf external-$TARGET/drm.tar.gz        && mv libdrm*           drm
    tar -xzf external-$TARGET/eigen.tar.gz      && mv eigen*          eigen
    tar -xzf external-$TARGET/expat.tar.gz      && mv expat*          expat
    tar -xzf external-$TARGET/libint.tar.gz     && mv libint*        libint
    tar -xzf external-$TARGET/mesa.tar.gz       && mv mesa*            mesa
    tar -xzf external-$TARGET/openblas.tar.gz   && mv OpenBLAS*    openblas
    tar -xzf external-$TARGET/pciaccess.tar.gz  && mv libpci*     pciaccess
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
    tar -xzf external-$TARGET/xshmfence.tar.gz  && mv libxshm*    xshmfence
    tar -xzf external-$TARGET/xtrans.tar.gz     && mv libxtrans*     xtrans
    tar -xzf external-$TARGET/xxf86vm.tar.gz    && mv libxxf86vm*   xxf86vm
    tar -xzf external-$TARGET/zlib.tar.gz       && mv zlib*            zlib

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
    cd xshmfence  && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..
    cd xxf86vm    && ./autogen.sh && PKG_CONFIG="pkg-config --static" ./configure "${CONFIGURE_GENERAL[@]}" && make -j $CORES && make install && cd ..

    cd expat && ./configure "${CONFIGURE_GENERAL[@]}"                        && make -j $CORES && make install && cd ..
    cd zlib  && ./configure --static --prefix="$DIRNAME/../external-$TARGET" && make -j $CORES && make install && cd ..

    cd pciaccess && meson setup build "${MESON_GENERAL[@]}"                   && meson compile -C build && meson install -C build && cd ..
    cd drm     && meson setup build "${MESON_GENERAL[@]}" "${MESON_DRM[@]}"   && meson compile -C build && meson install -C build && cd ..
    # cd mesa    && meson setup build "${MESON_GENERAL[@]}" "${MESON_MESA[@]}"  && meson compile -C build && meson install -C build && cd ..

    cd llvm   && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_OMP[@]}"    ../openmp && cmake_install && cd ../..
    cd eigen  && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}"                      ..        && cmake_install && cd ../..
    cd libint && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_LIBINT[@]}" ..        && cmake_install && cd ../..
    # cd raylib && mkdir -p build && cd build && cmake "${CMAKE_GENERAL[@]}" "${CMAKE_RAYLIB[@]}" ..        && cmake_install && cd ../..

    cd openblas && make "${MAKE_OPENBLAS[@]}" -j $CORES libs shared && make NO_SHARED=1 PREFIX="$DIRNAME/../external-$TARGET" install && cd ..

    rm -rf drm eigen expat libint llvm mesa openblas pciaccess raylib x11 xau xcb xcbproto xext xfixes xi xinerama xinput xorgmacros xorgproto xrender xrandr xshmfence xtrans xxf86vm zlib zigar zigcc zigcpp zigranlib

done
