LDhat Rust Bindings
===================

[English](./README.md)

## 依赖需求

1. [Rust 工具链 1.69.0](./rust-toolchain)
2. [Clang 5.0 或者更新](https://rust-lang.github.io/rust-bindgen/requirements.html) （可选）

推荐使用 [`rustup`](https://rustup.rs) 来设置 Rust 环境。

因为 LDhat 过时的 C 语法风格（大概是C89？），现代 Clang 例如 14.0 及其以后将无法编译 LDhat。
对源码进行简单的修改即可通过编译。或许可以更简单一点，直接注释掉[`Cargo.toml`](./Cargo.toml) 中的 `ldhat-sys` 依赖，因为目前项目并没有真正使用绑定，这样的话也就不需要 Clang 了。

https://github.com/tcztzy/ldhat-rs/blob/main/Cargo.toml#L17

## How to build

```console
$ git clone --recurse-submodules https://github.com/tcztzy/ldhat-rs
$ cd ldhat-rs
$ cargo build
```
