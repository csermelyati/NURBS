#!/bin/bash
name=(${1//\.*/})
g++ $1 -std=c++11 -lGL -lglut -lGLU -o $name
./$name
