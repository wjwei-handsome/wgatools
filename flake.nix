{
  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils, ... }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };
        rustPackage = pkgs.rustPlatform.buildRustPackage {
          pname = "wgatools";
          version = "0.0";
          src = self;
          cargoSha256 = "sha256-WvBJMBncDYaom54CJpxPU0805CKaYO8TXiSqxFAWKsg=";
          nativeBuildInputs = [ pkgs.cmake pkgs.pkg-config pkgs.openssl pkgs.perl ];
        };
        dockerImage = pkgs.dockerTools.buildImage {
          name = "wgatools";
          tag = "latest";
          copyToRoot = pkgs.buildEnv {
              name = "docker-root";
              paths = [ rustPackage ];
          };
        };
      in {
        packages = {
          wgatools = rustPackage;
          dockerImage = dockerImage;
        };
        defaultPackage = {
          x86_64-linux = rustPackage;
        };
        devShells.default = pkgs.mkShell {
          buildInputs = [ pkgs.rustc pkgs.cargo ];
        };
        dockerImage = {
          inherit dockerImage;
        };
      }
    );
}
