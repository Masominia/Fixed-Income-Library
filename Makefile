.PHONY: install build test lint run all

install:
	pip install -e ".[dev]"

CXX ?= g++
CXXFLAGS ?= -O2 -std=c++17 -Wall -Wextra -Icpp/include
LIB := $(filter-out cpp/src/main.cpp,$(wildcard cpp/src/*.cpp))

build:
	mkdir -p cpp/build
	$(CXX) $(CXXFLAGS) -o cpp/build/curvekit cpp/src/main.cpp $(LIB)
	$(CXX) $(CXXFLAGS) -o cpp/build/curvekit_tests $(wildcard cpp/tests/*.cpp) $(LIB)

lint:
	ruff check src tests scripts
	ruff format --check src tests scripts

test: build
	cpp/build/curvekit_tests
	pytest

run: build
	python -m curvekit run

all: lint test run
