"""
GROMACS Engine Interface
========================

Generate GROMACS simulation inputs for MDMix workflows.
Supports conversion from Amber topology and multi-chain systems.

Examples
--------
>>> from pymdmix.engines.gromacs import GromacsEngine
>>> engine = GromacsEngine()
>>> engine.write_minimization("system.top", "system.gro", "min.mdp")
"""
from __future__ import annotations

import logging
import re
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np

log = logging.getLogger(__name__)


@dataclass
class GromacsConfig:
    """
    GROMACS simulation configuration.
    
    Attributes
    ----------
    timestep : float
        Integration timestep in femtoseconds
    temperature : float
        Temperature in Kelvin
    pressure : float
        Pressure in bar
    cutoff : float
        Nonbonded cutoff in nm
    pme_spacing : float
        PME grid spacing in nm
    tau_t : float
        Temperature coupling time constant in ps
    tau_p : float
        Pressure coupling time constant in ps
    compressibility : float
        System compressibility in bar^-1
    num_threads : int
        Number of OpenMP threads
    gpu_id : str
        GPU device ID(s)
    """
    timestep: float = 2.0  # fs
    temperature: float = 300.0
    pressure: float = 1.0
    cutoff: float = 0.9  # nm
    pme_spacing: float = 0.12  # nm
    tau_t: float = 0.1  # ps
    tau_p: float = 2.0  # ps
    compressibility: float = 4.5e-5  # bar^-1
    num_threads: int = 0  # 0 = auto
    gpu_id: str = "0"
    
    # Restraint settings
    restraint_force: float = 1000.0  # kJ/mol/nm^2
    
    # Output frequencies (steps)
    nstlog: int = 5000
    nstxout: int = 0  # full precision coords
    nstvout: int = 0  # velocities
    nstfout: int = 0  # forces
    nstenergy: int = 5000
    nstxout_compressed: int = 500  # xtc output


class GromacsEngine:
    """
    GROMACS simulation engine for MDMix.
    
    Generates GROMACS MDP files and handles topology conversion
    from Amber format.
    """
    
    def __init__(self, config: GromacsConfig | None = None):
        self.config = config or GromacsConfig()
        self.log = logging.getLogger(self.__class__.__name__)
    
    def convert_amber_to_gromacs(
        self,
        prmtop: str | Path,
        inpcrd: str | Path,
        output_prefix: str = "system",
        output_dir: str | Path | None = None,
    ) -> tuple[Path, Path]:
        """
        Convert Amber topology to GROMACS format using parmed.
        
        Parameters
        ----------
        prmtop : Path
            Amber topology file
        inpcrd : Path
            Amber coordinate file
        output_prefix : str
            Prefix for output files
        output_dir : Path, optional
            Output directory
            
        Returns
        -------
        tuple[Path, Path]
            Paths to (topology.top, coordinates.gro)
        """
        prmtop = Path(prmtop)
        inpcrd = Path(inpcrd)
        output_dir = Path(output_dir) if output_dir else prmtop.parent
        
        try:
            import parmed as pmd
        except ImportError:
            raise ImportError("parmed required for Amber to GROMACS conversion")
        
        self.log.info(f"Converting {prmtop.name} to GROMACS format...")
        
        # Load Amber structure
        amber = pmd.load_file(str(prmtop), str(inpcrd))
        
        # Save GROMACS files
        top_path = output_dir / f"{output_prefix}.top"
        gro_path = output_dir / f"{output_prefix}.gro"
        
        amber.save(str(top_path), overwrite=True)
        amber.save(str(gro_path), overwrite=True)
        
        self.log.info(f"Created: {top_path.name}, {gro_path.name}")
        
        return top_path, gro_path
    
    def _get_common_mdp(self) -> str:
        """Generate common MDP parameters."""
        return f'''; Common parameters
integrator              = md
dt                      = {self.config.timestep / 1000:.4f}  ; ps
nstlog                  = {self.config.nstlog}
nstenergy               = {self.config.nstenergy}
nstxout-compressed      = {self.config.nstxout_compressed}
compressed-x-grps       = System

; Nonbonded
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rcoulomb                = {self.config.cutoff}
rvdw                    = {self.config.cutoff}

; Electrostatics
coulombtype             = PME
pme_order               = 4
fourierspacing          = {self.config.pme_spacing}

; Constraints
constraints             = h-bonds
constraint_algorithm    = lincs
lincs_iter              = 1
lincs_order             = 4

; Periodic boundary conditions
pbc                     = xyz
'''
    
    def _get_temperature_coupling(self, groups: str = "Protein Non-Protein") -> str:
        """Generate temperature coupling parameters."""
        group_list = groups.split()
        n_groups = len(group_list)
        temps = " ".join([str(self.config.temperature)] * n_groups)
        taus = " ".join([str(self.config.tau_t)] * n_groups)
        
        return f'''
; Temperature coupling
tcoupl                  = V-rescale
tc-grps                 = {groups}
tau_t                   = {taus}
ref_t                   = {temps}
'''
    
    def _get_pressure_coupling(self, coupling_type: str = "Parrinello-Rahman") -> str:
        """Generate pressure coupling parameters."""
        return f'''
; Pressure coupling
pcoupl                  = {coupling_type}
pcoupltype              = isotropic
tau_p                   = {self.config.tau_p}
ref_p                   = {self.config.pressure}
compressibility         = {self.config.compressibility}
'''
    
    def write_minimization(
        self,
        output: str | Path,
        steps: int = 5000,
        algorithm: str = "steep",
        restraints: bool = True,
    ) -> Path:
        """
        Write GROMACS minimization MDP file.
        
        Parameters
        ----------
        output : Path
            Output MDP file path
        steps : int
            Maximum minimization steps
        algorithm : str
            Minimization algorithm (steep, cg, l-bfgs)
        restraints : bool
            Include position restraints
            
        Returns
        -------
        Path
            Path to generated MDP file
        """
        output = Path(output)
        
        define = "-DPOSRES" if restraints else ""
        
        mdp = f'''; Minimization parameters
define                  = {define}
integrator              = {algorithm}
nsteps                  = {steps}
emtol                   = 100.0
emstep                  = 0.01

; Nonbonded
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 10
rcoulomb                = {self.config.cutoff}
rvdw                    = {self.config.cutoff}

; Electrostatics
coulombtype             = PME
pme_order               = 4
fourierspacing          = {self.config.pme_spacing}

; Output
nstxout                 = 0
nstvout                 = 0
nstenergy               = 500
nstlog                  = 500

; PBC
pbc                     = xyz
'''
        
        output.write_text(mdp)
        self.log.info(f"Wrote minimization MDP: {output}")
        return output
    
    def write_nvt_equilibration(
        self,
        output: str | Path,
        steps: int = 50000,
        restraints: bool = True,
        gen_vel: bool = True,
        temperature: float | None = None,
        tc_groups: str = "Protein Non-Protein",
    ) -> Path:
        """
        Write GROMACS NVT equilibration MDP file.
        
        Parameters
        ----------
        output : Path
            Output MDP file path
        steps : int
            Number of equilibration steps
        restraints : bool
            Include position restraints
        gen_vel : bool
            Generate initial velocities
        temperature : float, optional
            Temperature (defaults to config value)
        tc_groups : str
            Temperature coupling groups
            
        Returns
        -------
        Path
            Path to generated MDP file
        """
        output = Path(output)
        temperature = temperature or self.config.temperature
        
        define = "-DPOSRES" if restraints else ""
        gen_vel_str = "yes" if gen_vel else "no"
        
        mdp = f'''; NVT Equilibration
define                  = {define}
{self._get_common_mdp()}
nsteps                  = {steps}

; Velocity generation
gen_vel                 = {gen_vel_str}
gen_temp                = {temperature}
gen_seed                = -1

{self._get_temperature_coupling(tc_groups)}

; No pressure coupling in NVT
pcoupl                  = no
'''
        
        output.write_text(mdp)
        self.log.info(f"Wrote NVT equilibration MDP: {output}")
        return output
    
    def write_npt_equilibration(
        self,
        output: str | Path,
        steps: int = 100000,
        restraints: bool = True,
        tc_groups: str = "Protein Non-Protein",
    ) -> Path:
        """
        Write GROMACS NPT equilibration MDP file.
        
        Parameters
        ----------
        output : Path
            Output MDP file path
        steps : int
            Number of equilibration steps
        restraints : bool
            Include position restraints
        tc_groups : str
            Temperature coupling groups
            
        Returns
        -------
        Path
            Path to generated MDP file
        """
        output = Path(output)
        
        define = "-DPOSRES" if restraints else ""
        
        mdp = f'''; NPT Equilibration
define                  = {define}
{self._get_common_mdp()}
nsteps                  = {steps}

; Continue from NVT
gen_vel                 = no
continuation            = yes

{self._get_temperature_coupling(tc_groups)}
{self._get_pressure_coupling("Berendsen")}
refcoord_scaling        = com
'''
        
        output.write_text(mdp)
        self.log.info(f"Wrote NPT equilibration MDP: {output}")
        return output
    
    def write_production(
        self,
        output: str | Path,
        steps: int = 500000,
        ensemble: str = "NVT",
        tc_groups: str = "Protein Non-Protein",
    ) -> Path:
        """
        Write GROMACS production MD MDP file.
        
        Parameters
        ----------
        output : Path
            Output MDP file path
        steps : int
            Number of production steps
        ensemble : str
            Ensemble type (NPT or NVT)
        tc_groups : str
            Temperature coupling groups
            
        Returns
        -------
        Path
            Path to generated MDP file
        """
        output = Path(output)
        
        pressure_section = ""
        if ensemble == "NPT":
            pressure_section = self._get_pressure_coupling("Parrinello-Rahman")
        else:
            pressure_section = "\n; No pressure coupling (NVT)\npcoupl                  = no\n"
        
        mdp = f'''; Production MD ({ensemble})
{self._get_common_mdp()}
nsteps                  = {steps}

; Continue from equilibration
gen_vel                 = no
continuation            = yes

{self._get_temperature_coupling(tc_groups)}
{pressure_section}
'''
        
        output.write_text(mdp)
        self.log.info(f"Wrote production MDP ({ensemble}): {output}")
        return output
    
    def write_restraints_itp(
        self,
        gro_file: str | Path,
        output: str | Path,
        selection: str = "protein",
        force_constant: float | None = None,
    ) -> Path:
        """
        Write position restraints ITP file.
        
        Parameters
        ----------
        gro_file : Path
            GRO coordinate file to extract atom indices
        output : Path
            Output ITP file path
        selection : str
            Atom selection (protein, backbone, heavy)
        force_constant : float, optional
            Force constant in kJ/mol/nm^2
            
        Returns
        -------
        Path
            Path to generated ITP file
        """
        gro_file = Path(gro_file)
        output = Path(output)
        force_constant = force_constant or self.config.restraint_force
        
        # Parse GRO file to get atom indices
        atoms = []
        with open(gro_file) as f:
            lines = f.readlines()
            n_atoms = int(lines[1].strip())
            for i, line in enumerate(lines[2:2+n_atoms], start=1):
                resname = line[5:10].strip()
                atomname = line[10:15].strip()
                
                if selection == "protein":
                    # Include all protein heavy atoms
                    if resname not in ["WAT", "HOH", "SOL", "NA", "CL", "Na+", "Cl-"]:
                        if not atomname.startswith("H"):
                            atoms.append(i)
                elif selection == "backbone":
                    if atomname in ["CA", "C", "N", "O"]:
                        atoms.append(i)
                elif selection == "heavy":
                    if not atomname.startswith("H"):
                        atoms.append(i)
        
        # Write ITP
        fc = int(force_constant)
        lines = [
            "; Position restraints",
            "[ position_restraints ]",
            "; atom  type      fx      fy      fz",
        ]
        for atom_idx in atoms:
            lines.append(f"   {atom_idx:5d}     1  {fc}  {fc}  {fc}")
        
        output.write_text("\n".join(lines))
        self.log.info(f"Wrote restraints ITP: {output} ({len(atoms)} atoms)")
        return output
    
    def create_index_groups(
        self,
        gro_file: str | Path,
        output: str | Path,
        solvent_residues: list[str] | None = None,
    ) -> Path:
        """
        Create GROMACS index file with custom groups.
        
        Parameters
        ----------
        gro_file : Path
            Input GRO file
        output : Path
            Output NDX file
        solvent_residues : list, optional
            Solvent residue names (e.g., ["WAT", "ETA"])
            
        Returns
        -------
        Path
            Path to generated NDX file
        """
        gro_file = Path(gro_file)
        output = Path(output)
        solvent_residues = solvent_residues or ["WAT", "SOL", "HOH"]
        
        # Use gmx make_ndx if available
        try:
            cmd = f'echo "q" | gmx make_ndx -f {gro_file} -o {output}'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode == 0:
                self.log.info(f"Created index file: {output}")
                return output
        except Exception:
            pass
        
        # Fallback: create basic index file manually
        self.log.warning("gmx make_ndx not available, creating basic index")
        
        # Parse GRO and create groups
        protein_atoms = []
        solvent_atoms = []
        system_atoms = []
        
        with open(gro_file) as f:
            lines = f.readlines()
            n_atoms = int(lines[1].strip())
            for i, line in enumerate(lines[2:2+n_atoms], start=1):
                resname = line[5:10].strip()
                system_atoms.append(i)
                
                if resname in solvent_residues or resname in ["NA", "CL", "Na+", "Cl-"]:
                    solvent_atoms.append(i)
                else:
                    protein_atoms.append(i)
        
        # Write NDX
        with open(output, "w") as f:
            f.write("[ System ]\n")
            f.write(" ".join(map(str, system_atoms)) + "\n")
            f.write("[ Protein ]\n")
            f.write(" ".join(map(str, protein_atoms)) + "\n")
            f.write("[ Non-Protein ]\n")
            f.write(" ".join(map(str, solvent_atoms)) + "\n")
        
        self.log.info(f"Created index file: {output}")
        return output
    
    def prealign_trajectory(
        self,
        tpr_file: str | Path,
        xtc_file: str | Path,
        output: str | Path | None = None,
    ) -> Path:
        """
        Pre-align trajectory for analysis (center and image).
        
        Uses gmx trjconv to center on protein and apply PBC corrections.
        
        Parameters
        ----------
        tpr_file : Path
            TPR run file
        xtc_file : Path
            Input XTC trajectory
        output : Path, optional
            Output XTC (defaults to input with _aligned suffix)
            
        Returns
        -------
        Path
            Path to aligned trajectory
        """
        tpr_file = Path(tpr_file)
        xtc_file = Path(xtc_file)
        output = Path(output) if output else xtc_file.with_stem(xtc_file.stem + "_aligned")
        
        tmp_file = xtc_file.with_stem(xtc_file.stem + "_tmp")
        
        # Step 1: Remove jumps
        cmd1 = f"echo '1 0' | gmx trjconv -s {tpr_file} -f {xtc_file} -o {tmp_file} -pbc nojump -center"
        
        # Step 2: Apply compact representation
        cmd2 = f"echo '1 0' | gmx trjconv -s {tpr_file} -f {tmp_file} -o {output} -pbc mol -ur compact -center"
        
        self.log.info(f"Pre-aligning trajectory: {xtc_file.name}")
        
        try:
            subprocess.run(cmd1, shell=True, check=True, capture_output=True)
            subprocess.run(cmd2, shell=True, check=True, capture_output=True)
            tmp_file.unlink(missing_ok=True)
            self.log.info(f"Aligned trajectory: {output}")
            return output
        except subprocess.CalledProcessError as e:
            self.log.error(f"Trajectory alignment failed: {e}")
            raise
    
    def write_replica_inputs(
        self,
        replica,
        output_dir: str | Path | None = None,
        convert_topology: bool = True,
    ) -> list[Path]:
        """
        Write all GROMACS input files for a replica.
        
        Parameters
        ----------
        replica : Replica
            Replica to generate inputs for
        output_dir : Path, optional
            Output directory (defaults to replica.path)
        convert_topology : bool
            Convert Amber topology to GROMACS format
            
        Returns
        -------
        list[Path]
            List of generated file paths
        """
        output_dir = Path(output_dir) if output_dir else replica.path
        output_dir.mkdir(parents=True, exist_ok=True)
        
        files = []
        
        # Convert topology if needed
        if convert_topology and replica.topology.endswith(".prmtop"):
            top_path, gro_path = self.convert_amber_to_gromacs(
                replica.path / replica.topology,
                replica.path / replica.coordinates,
                output_prefix="system",
                output_dir=output_dir,
            )
            files.extend([top_path, gro_path])
        else:
            gro_path = output_dir / "system.gro"
        
        # Create index groups
        ndx_path = self.create_index_groups(gro_path, output_dir / "index.ndx")
        files.append(ndx_path)
        
        # Create restraints
        posre_path = self.write_restraints_itp(gro_path, output_dir / "posre.itp")
        files.append(posre_path)
        
        # Minimization (2 stages)
        files.append(self.write_minimization(output_dir / "min1.mdp", steps=5000, algorithm="steep"))
        files.append(self.write_minimization(output_dir / "min2.mdp", steps=5000, algorithm="cg"))
        
        # NVT equilibration
        files.append(self.write_nvt_equilibration(output_dir / "nvt.mdp", steps=50000))
        
        # NPT equilibration
        files.append(self.write_npt_equilibration(output_dir / "npt.mdp", steps=100000))
        
        # Production
        files.append(self.write_production(output_dir / "prod.mdp", steps=500000, ensemble="NVT"))
        
        # Write run script
        run_script = output_dir / "run_gromacs.sh"
        run_script.write_text(f'''#!/bin/bash
# GROMACS run script - Generated by pyMDMix

GMX="gmx"
TOP="system.top"
GRO="system.gro"
NDX="index.ndx"

# Minimization 1 (steepest descent)
$GMX grompp -f min1.mdp -c $GRO -p $TOP -o min1.tpr -maxwarn 2
$GMX mdrun -deffnm min1 -v

# Minimization 2 (conjugate gradient)
$GMX grompp -f min2.mdp -c min1.gro -p $TOP -o min2.tpr -maxwarn 2
$GMX mdrun -deffnm min2 -v

# NVT equilibration
$GMX grompp -f nvt.mdp -c min2.gro -r min2.gro -p $TOP -n $NDX -o nvt.tpr -maxwarn 2
$GMX mdrun -deffnm nvt -v

# NPT equilibration
$GMX grompp -f npt.mdp -c nvt.gro -r nvt.gro -p $TOP -n $NDX -o npt.tpr -maxwarn 2
$GMX mdrun -deffnm npt -v

# Production
$GMX grompp -f prod.mdp -c npt.gro -p $TOP -n $NDX -o prod.tpr -maxwarn 2
$GMX mdrun -deffnm prod -v

echo "Done!"
''')
        run_script.chmod(0o755)
        files.append(run_script)
        
        return files


@dataclass
class GromacsCheckResult:
    """Result of GROMACS simulation check."""
    minimization: bool = False
    equilibration: bool = False
    production: bool = False
    
    @property
    def complete(self) -> bool:
        return self.minimization and self.equilibration and self.production


class GromacsCheck:
    """
    Check GROMACS simulation status.
    
    Verifies that simulation stages completed successfully.
    """
    
    def __init__(self):
        self.log = logging.getLogger(self.__class__.__name__)
    
    def check_minimization(self, path: str | Path) -> bool:
        """Check if minimization completed."""
        path = Path(path)
        log_file = path / "min2.log"
        
        if not log_file.exists():
            log_file = path / "min1.log"
        
        if not log_file.exists():
            return False
        
        content = log_file.read_text()
        return "Finished mdrun" in content or "Force converged" in content
    
    def check_equilibration(self, path: str | Path) -> bool:
        """Check if equilibration completed."""
        path = Path(path)
        
        for stage in ["nvt", "npt"]:
            log_file = path / f"{stage}.log"
            if not log_file.exists():
                return False
            
            content = log_file.read_text()
            if "Finished mdrun" not in content:
                return False
        
        return True
    
    def check_production(self, path: str | Path) -> bool:
        """Check if production completed."""
        path = Path(path)
        log_file = path / "prod.log"
        
        if not log_file.exists():
            return False
        
        content = log_file.read_text()
        return "Finished mdrun" in content
    
    def check_all(self, path: str | Path) -> GromacsCheckResult:
        """Check all simulation stages."""
        return GromacsCheckResult(
            minimization=self.check_minimization(path),
            equilibration=self.check_equilibration(path),
            production=self.check_production(path),
        )
    
    def get_box_volume(self, log_file: str | Path) -> float | None:
        """
        Extract average box volume from production log.
        
        Parameters
        ----------
        log_file : Path
            GROMACS log file
            
        Returns
        -------
        float or None
            Box volume in Å³, or None if not found
        """
        log_file = Path(log_file)
        
        if not log_file.exists():
            return None
        
        content = log_file.read_text()
        
        # Look for box dimensions in log
        # Box-X, Box-Y, Box-Z in nm
        match = re.search(r"Box-X\s+Box-Y\s+Box-Z\s*\n\s*([\d.e+-]+)\s+([\d.e+-]+)\s+([\d.e+-]+)", content)
        
        if match:
            box_x = float(match.group(1)) * 10  # nm to Å
            box_y = float(match.group(2)) * 10
            box_z = float(match.group(3)) * 10
            # Apply truncated octahedron correction if applicable
            volume = box_x * box_y * box_z * 0.77
            return volume
        
        return None
