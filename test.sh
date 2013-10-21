#!/bin/sh
export RUST_TEST_TASKS=1 
rustc --test wgfem.rs && ./wgfem
