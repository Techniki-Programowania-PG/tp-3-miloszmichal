cmake_minimum_required(VERSION 3.19)

project(scikit_build_example LANGUAGES CXX)



if (MSVC)
  set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} /WX-")
  set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} /WX-")
  set(CMAKE_CXX_FLAGS_RELEASE"${CMAKE_CXX_FLAGS_RELEASE} /WX-")
else()
  set(CMAKE_CXX_FLAGS        "${CMAKE_CXX_FLAGS} -Wno-error")
  set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -Wno-error")
  set(CMAKE_CXX_FLAGS_RELEASE"${CMAKE_CXX_FLAGS_RELEASE} -Wno-error")
endif()

# 
# FetchContent: pobieramy pybind11 i matplotplusplus
# 
include(FetchContent)

# pybind11 (bindingi C++ ↔ Python)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.12.0
)
FetchContent_MakeAvailable(pybind11)

#  matplotplusplus (wizualizacja wykresów)
FetchContent_Declare(
  matplotplusplus
  GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus.git
  GIT_TAG        origin/master
)
FetchContent_MakeAvailable(matplotplusplus)

# 
# Usuwamy /WX i -Werror z targetu matplot, jeśli Matplot++ go dodał
# 
if (TARGET matplot)
  if (MSVC)
    # Dodanie "/WX-" do kompilacji matplot, aby na wszelki wypadek wyłączyć
    # traktowanie ostrzeżeń jako błędów (jeśli Matplot++ dodał /WX).
    target_compile_options(matplot PUBLIC "/WX-")
  else()
    # Dodanie "-Wno-error", aby wyłączyć -Werror w matplot
    target_compile_options(matplot PUBLIC "-Wno-error")
  endif()
endif()

# 
# Biblioteka DSP (statyczna) z plików w src/dsp/
# 
add_library(dsp STATIC
    src/dsp/signal.cpp
)

# Ponieważ w bindingach używamy `#include "signal.hpp"`, dodajemy do include path
# bezpośrednio katalog src/dsp, aby kompilator znalazł `signal.hpp`.
target_include_directories(dsp PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp
)

# Matplotplusplus daje target `matplot`, z którego łączymy `dsp`.
target_link_libraries(dsp PUBLIC matplot)

#
# Moduł Python scikit_build_example (pybind11)
# 
pybind11_add_module(scikit_build_example
    src/scikit_build_example/module.cpp
    src/scikit_build_example/generators.cpp
    src/scikit_build_example/transforms.cpp
    src/scikit_build_example/filters.cpp
    src/scikit_build_example/plot.cpp
    src/scikit_build_example/utils.cpp 
)

# Linkujemy statyczną bibliotekę dsp do modułu Python
target_link_libraries(scikit_build_example PRIVATE dsp)

# 
# Wymuszenie standardu C++17 (Matplotplusplus wymaga co najmniej C++17)
set_target_properties(dsp PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)
set_target_properties(scikit_build_example PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
)
