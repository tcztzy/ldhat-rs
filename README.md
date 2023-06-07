LDhat Rust Bindings
===================

[中文](./README_zh_CN.md)

## Requirement

1. [Rust toolchain 1.69.0](./rust-toolchain)
2. [Clang 5.0 or later](https://rust-lang.github.io/rust-bindgen/requirements.html) (Optional)

It is recommend use [`rustup`](https://rustup.rs) to setup your Rust environment.

TIPS:

modern clang such as 14.0 or later will not compile the LDhat because the outdated C code style (maybe C89?). It is very easy manually patch the source codes to make it compiled. Or even simpler, just comment out the `ldhat-sys` dependency in [`Cargo.toml`](./Cargo.toml) since the project did not use the binding actually.

https://github.com/tcztzy/ldhat-rs/blob/main/Cargo.toml#L17

## How to build

```console
$ git clone --recurse-submodules https://github.com/tcztzy/ldhat-rs
$ cd ldhat-rs
$ cargo build
```
