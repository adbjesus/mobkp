{
  inputs = {
    nixpkgs = {
      url = "github:nixos/nixpkgs/nixos-unstable";
    };

    flake-utils = {
      url = "github:numtide/flake-utils";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    mooutils = {
      url = "github:adbjesus/mooutils";
      inputs.nixpkgs.follows = "nixpkgs";
    };

    apm = {
      url = "github:adbjesus/apm";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs = { self, nixpkgs, flake-utils, mooutils, apm }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = nixpkgs.legacyPackages.${system};
      in rec {
        packages.mobkp = pkgs.stdenv.mkDerivation {
          pname = "mobkp";
          version = "0.1.0";
          src = self;

          meta = with nixpkgs.lib; {
            description = "Algorithms for the Multi-Objective Knapsack Problem (MOBKP)";
            license = licenses.mit;
          };

          nativeBuildInputs = with pkgs; [
            cmake
            ninja
          ];

          buildInputs = [
            pkgs.glpk
            pkgs.fmt_8
            pkgs.boost175
            pkgs.cli11
            mooutils.packages.${system}.mooutils
            apm.packages.${system}.apm
          ];
        };

        packages.default = self.packages.${system}.mobkp;
      });
}
