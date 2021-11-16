{
  description = "mobkp";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs";
    flake-utils.url = "github:numtide/flake-utils";
    mooutils.url = "git+ssh://git@git.adbjesus.com/mooutils";
    apm.url = "git+ssh://git@github.com/adbjesus/apm";
  };

  outputs = { self, nixpkgs, flake-utils, mooutils, apm }:
    flake-utils.lib.eachDefaultSystem (system:
      let pkgs = nixpkgs.legacyPackages.${system};
      in rec {
        packages.mobkp = pkgs.stdenv.mkDerivation {
          pname = "mobkp";
          version = "0.1.0";
          src = self;
          nativeBuildInputs = with pkgs; [ cmake ninja ];
          buildInputs = [
            pkgs.glpk
            pkgs.fmt_8
            pkgs.catch2
            pkgs.cli11
            pkgs.boost175
            mooutils.packages.${system}.mooutils
            apm.packages.${system}.apm
          ];
        };
        defaultPackage = self.packages.${system}.mobkp;
      });
}
