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
          name = "wgatools";
          src = self;
          cargoSha256 = "0000000000000000000000000000000000000000000000000000";
        };
        dockerImage = pkgs.dockerTools.buildImage {
          name = "wgatools";
          tag = "latest";
          contents = [ rustPackage ];
        };
      in {
        packages = {
          inherit rustPackage;
        };
        defaultPackage = {
          inherit rustPackage;
        };
        devShells.default = pkgs.mkShell {
          buildInputs = [ pkgs.rustc pkgs.cargo ];
        };
        dockerImage = {
          inherit dockerImage;
        };
        // Add the Singularity image here
        singularityImage = {
          inherit singularityImage;
        };
      }
    );
}
