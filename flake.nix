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
          version = "1.0";  // Replace with your package's version
          src = self;
          cargoSha256 = "0000000000000000000000000000000000000000000000000000"; // Replace after the first build attempt
        };
        dockerImage = pkgs.dockerTools.buildImage {
          name = "wgatools";
          tag = "latest";
          contents = [ rustPackage ];
        };
      in {
        packages = {
          wgatools = rustPackage;
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
