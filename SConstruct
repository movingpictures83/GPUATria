#!/usr/bin/env python3
# Copyright (C) 2021 FIUBioRG
# SPDX-License-Identifier: MIT
#
# -*-python-*-

import logging
import os
from os import environ, getenv
from os.path import abspath, realpath, dirname
import sys
import subprocess

if sys.version_info.major < 3:
    logging.error("Python3 is required to run this SConstruct")

AddOption(
    "--pluma",
    dest="pluma",
    type="string",
    nargs=1,
    action="store",
    metavar="DIR",
    default=abspath("../../src"),
    help="DIR for the installation of PluMA's header files"
)

AddOption(
    "--gpu-architecture",
    dest="gpu_arch",
    type="string",
    nargs=1,
    action="store",
    metavar="ARCH",
    default="sm_35",
    help="Specify the name of the class of NVIDIA 'virtual' GPU architecture for which the CUDA input files must be compiled"
)

env = Environment(
    Env=environ,
    PATH=getenv("path"),
    CC=getenv("CC", "cc"),
    CXX=getenv("CXX", "c++"),
    NVCC=getenv("NVCC", "nvcc"),
    SHCCFLAGS=["-fpermissive", "-fPIC", "-I.", "-O2"],
    SHCXXFLAGS=["-std=c++11", "-fpermissive", "-fPIC", "-I.", "-O2"],
    CCFLAGS=["-fpermissive", "-fPIC", "-I.", "-O2"],
    CXXFLAGS=["-std=c++11", "-fpermissive", "-fPIC", "-I.", "-O2"],
    LICENSE=["MIT"],
    SHLIBPREFIX="lib",
    CUDA_PATH=getenv("CUDA_PATH", "/usr/local/cuda"),
    CUDA_SDK_PATH=getenv("CUDA_SDK_PATH", "/usr/local/cuda"),
    NVCCFLAGS=[
      "-I" + os.getcwd(),
      "--ptxas-options=-v",
      "-std=c++11",
      "-Xcompiler",
      "-fPIC"
    ],
    GPU_ARCH=GetOption("gpu_arch")
)

env.Append(
    SHCXXFLAGS=[
      "-I" + env.GetOption("pluma"),
      "-I" + str(env['CUDA_PATH']) + '/include'
    ],
    CXXFLAGS=[
      "-I" + env.GetOption("pluma"),
      "-I" + str(env['CUDA_PATH']) + '/include'
    ],
    NVCCFLAGS=[
      "-I" + env.GetOption("pluma"),
      "-I" + str(env['CUDA_PATH']) + '/include'
    ]
)

if not sys.platform.startswith("darwin"):
    env.Append(LINKFLAGS=["-rdynamic"])
    env.Append(LIBS=["rt"])
else:
    env.Append(CCFLAGS=["-DAPPLE"])

config = Configure(env)

if not config.CheckCC():
    Exit(1)

if not config.CheckCXX():
    Exit(1)

if not config.CheckSHCXX():
    Exit(1)

if not config.CheckCXXHeader("csv_parser/csv_parser.hpp"):
    Exit(1)

libs = [
    "m",
    "c"
]

for lib in libs:
    if not config.CheckLib(lib):
        Exit(1)

plugin = realpath("./GPUAtriaPlugin.cu")
filesInPath = Glob(dirname(plugin) + "./*.cu")
pluginName = plugin[plugin.rindex("/") + 1:]
pluginName = pluginName.replace(".cu", "")
pluginName = "lib" + pluginName

env.Command(
  pluginName,
  filesInPath,
  "nvcc -o $TARGET -shared $NVCCFLAGS -arch=$GPU_ARCH $SOURCES"
)