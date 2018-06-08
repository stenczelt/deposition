program prep_for_impact

! This program takes (current) output LAMMPS structure files and
! generates an input file with a randomly (XY) placed incident ion at a suitable (Z) distance
! from the substrate

  use neighbors

  implicit none

  integer :: i, n_atoms, i_crap, j, k, n_remove, n_clusters, largest_cluster
  character*64 :: output, c_crap, input
  real*8 :: lx(1:2), ly(1:2), lz(1:2), ion_energy, pos(1:3), vel(1:3), dz, dxy, rand
  real*8, allocatable :: positions(:,:), velocities(:,:)
  integer :: n, clock, cluster
  integer, allocatable :: seed(:), atom_belongs_to_cluster(:), cluster_size(:)
  logical :: add_ion = .true.
  real*8 :: padding, d, bonding_cutoff
  logical, allocatable :: bonded(:,:), remove_atom(:), clustered(:)
  logical, allocatable :: atom_visited(:)
  integer :: iostatus, nlines
  character*128 :: line, option, value
  character*1 :: comment_line
  real*8 :: ion_mass, min_dist




! Input parameters and options
!
! Beware of nomenclature. "output" stands for output of the previous LAMMPS simulation
! and "input" stands for input for the subsequent LAMMPS simulation. *This* code takes
! "output" as its input and "input" is its output
! Defaults
  output = "output_structure"
  input = "input_structure"
  ion_energy = 0.d0
  bonding_cutoff = 1.9d0
  padding = 4.d0
  ion_mass = 12.01d0
  min_dist = 3.d0
! If the user gives an "options" file, we override the defaults accordingly
  open(unit=10, file="prep_for_impact_options", status="old")
  iostatus = 0
! NOTE: this input reader is not fool proof and will break easily. Make sure you don't
! expect too much from the code when you use it
  do while( iostatus == 0 )
    read(10,"(A)",iostat=iostatus) line
    comment_line = line
    if( comment_line /= "!" .and. comment_line /= "#" )then
      read(line, *) option, value, value
      if( trim(adjustl(option)) == "output" )then
        read(value, *) output
      else if( trim(adjustl(option)) == "input" )then
        read(value, *) input
      else if( trim(adjustl(option)) == "ion_energy" )then
        read(value, *) ion_energy
      else if( trim(adjustl(option)) == "bonding_cutoff" )then
        read(value, *) bonding_cutoff
      else if( trim(adjustl(option)) == "padding" )then
        read(value, *) padding
      else if( trim(adjustl(option)) == "ion_mass" )then
        read(value, *) ion_mass
      else if( trim(adjustl(option)) == "min_dist" )then
        read(value, *) min_dist
      end if
    end if
  end do



  write(*,*) "____________________________________________________________________"
  write(*,*) "                                                                    \"
  write(*,*) "Welcome to Miguel Caro's 'preparation for impact' code.             |"
  write(*,*) "This code takes LAMMPS output trajectory files as input, it adds an |"
  write(*,*) "impacting atom and prepares the structure for subsequent LAMMPS     |"
  write(*,*) "MD simulation. The code is meant for deposition simulations similar |"
  write(*,*) "to those reported in Phys. Rev. Lett. 120, 166101 (2018).           |"
  write(*,*) "                                                                    |"
  write(*,*) "The code is copyright of Miguel Caro (mcaroba@gmail.com) and was    |"
  write(*,*) "written during his employment at Aalto University in the period     |"
  write(*,*) "2017-2018. For enquiries you can contact the author by email.       |"
  write(*,*) "                                                                    |"
  write(*,*) "This code was last modified on 08 Jun 2018                          |"
  write(*,*) "                                                    ________________/"
  write(*,*) "                                                   /"
  write(*,*) "                                                   |"
  write(*,*) "***************************************************|"
  write(*,*) "                                                   |"
  write(*,*) "You have chosen the following options:             |"
  write(*,*) "                                                   |"
  write(*,*) "*LAMMPS output file from previous simulation:      |"
  if( len(trim(adjustl(output))) > 49 )then
    write(*,"(A, A44, A)")  "    ", trim(adjustl(output)), "... |"
  else
    write(*,"(A, A47, A)")  "    ", trim(adjustl(output)), " |"
  end if
  write(*,*) "                                                   |"
  write(*,*) "*Input file for subsequent LAMMPS simulation:      |"
  if( len(trim(adjustl(input))) > 49 )then
    write(*,"(A, A44, A)")  "    ", trim(adjustl(input)), "... |"
  else
    write(*,"(A, A47, A)")  "    ", trim(adjustl(input)), " |"
  end if
  write(*,*) "                                                   |"
  write(*,"(A, F24.3, A)") " *Incident atom energy: ", ion_energy, " eV |"
  write(*,*) "                                                   |"
  write(*,"(A, F25.3, A)") " *Incident atom mass: ", ion_mass, " amu |"
  write(*,*) "                                                   |"
  write(*,"(A, F19.3, A)") " *Bonded cutoff (1st NN):", bonding_cutoff, " Angst. |"
  write(*,*) "                                                   |"
  write(*,"(A, F20.3, A)") " *Cell padding (vacuum):", padding, " Angst. |"
  write(*,*) "                                                   |"


! If ion energy is zero, then we just transform the LAMMPS output to LAMMPS input
! and do not add any incident atoms
  if( ion_energy == 0.d0 )then
    add_ion=.false.
  end if


! Read previous LAMMPS output file
! We only need to read the *last* snapshot of the trajectory. This code does not assume
! any particular number of snapshots, or a constant atom number, for that matter. It will
! detect which is the last snapshot and how many atoms are contained in it, and will use
! that configuration for the processing.
  nlines = 0
  open(unit=10, file=output, status="old")
  read(10, "(A)", iostat=iostatus) line
  if( trim(adjustl(line)) /= "ITEM: TIMESTEP" )then
    write(*,*) "***************************************************|"
    write(*,*) "                                                   |"
    write(*,*) "ERROR:                                             | <-- ERROR"
    write(*,*) "                                                   |"
    write(*,*) "The file used does not seem to be a native LAMMPS  |"
    write(*,*) "output trajectory file.                            |"
    write(*,*) "                                                   |"
    write(*,*) "___________________________________________________/"
    return
  end if
  nlines = nlines + 1
  do while( iostatus == 0 )
    read(10,"(A)", iostat=iostatus) line
    nlines = nlines + 1
    if( trim(adjustl(line)) == "ITEM: NUMBER OF ATOMS" )then
      read(10,*, iostat=iostatus) n_atoms
      do i = 1, n_atoms + 5
        read(10, *, iostat=iostatus)
      end do
      nlines = nlines + n_atoms + 6
      read(10,"(A)", iostat=iostatus)
      backspace(10)
    end if
  end do
  allocate( positions(1:n_atoms, 1:3) )
  allocate( velocities(1:n_atoms, 1:3) )
  rewind(10)
  write(*,*) "***************************************************|"
  write(*,*) "                                                   |"
  write(*,"(A, I8, A)") " Reading positions/velocities for ", n_atoms, " atoms... |"
  write(*,*) "                                                   |"
  do i = 1, nlines - n_atoms - 4
    read(10, *)
  end do
  read(10,*) lx(1:2)
  read(10,*) ly(1:2)
  read(10,*) lz(1:2)
  read(10,*)
  do i = 1, n_atoms
    read(10,*) c_crap, positions(i, 1:3), velocities(i, 1:3)
  end do
  close(10)




! Make sure that atoms are not drifting away from the slab. I have observed this
! for instance with the Tersoff potential. Then two-atom clusters start filling
! the vacuum region and mess up the simulation. What we do here is identify atoms
! or clusters of atoms which are not connected to the rest of
! the slab and we remove them from the simulation box
  write(*,*) "***************************************************|"
  write(*,*) "                                                   |"
  write(*,*) "Making sure that chunks of the growing film are    |"
  write(*,*) "not floating away... (floating chunks will be      |"
  write(*,*) "removed from the simulation box)                   |"
  write(*,*) "                                                   |"
  allocate( bonded(1:n_atoms, 1:n_atoms) )
  allocate( remove_atom(1:n_atoms) )
  bonded = .false.
  remove_atom = .false.
! Check if two atoms are bonded and whether an atom is clustered
!$omp parallel do private(i,j,d)
  do i = 1, n_atoms
    do j = i+1, n_atoms
      call get_distance(positions(i,1:3), positions(j,1:3), lx, ly, lz, d)
      if( d < bonding_cutoff )then
        bonded(i, j) = .true.
        bonded(j, i) = .true.
      end if
    end do
  end do
! Check whether an atom is clustered.
! We build clusters of atoms depending on whether we can connect an ensemble of atoms
! through bonds. There will be at least one big cluster (the growing film). If there
! are other (isolated) clusters besides the larger one, then we remove them
  allocate( atom_belongs_to_cluster(1:n_atoms) )
  atom_belongs_to_cluster = 0
  allocate( atom_visited(1:n_atoms) )
  atom_visited = .false.
  cluster = 0
  do i = 1, n_atoms
!   Find all the atoms that can be connected to i and put them on the same cluster
    if( .not. atom_visited(i) )then
      cluster = cluster + 1
      atom_belongs_to_cluster(i) = cluster
      atom_visited(i) = .true.
      do j = 1, n_atoms
        call find_neighbors(i, j, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)
      end do
    end if
  end do
  n_clusters = cluster
! Calculate cluster sizes
  allocate( cluster_size(1:n_clusters) )
  cluster_size = 0
  do i = 1, n_atoms
    j = atom_belongs_to_cluster(i)
    cluster_size(j) = cluster_size(j) + 1
  end do
  largest_cluster = 1
  do i = 1, n_clusters
    if( cluster_size(i) > cluster_size(largest_cluster) )then
      largest_cluster = i
    end if
  end do
! Mark atoms which are not in the main cluster for deletion
  do i = 1, n_atoms
    if( atom_belongs_to_cluster(i) /= largest_cluster )then
      remove_atom(i) = .true.
    end if
  end do
! Check number of atoms to be removed
  n_remove = 0
  do i = 1, n_atoms
    if( remove_atom(i) )then
      n_remove = n_remove + 1
    end if
  end do





  write(*,*) "***************************************************|"
  write(*,*) "                                                   |"
  write(*,*) "Randomizing XY position of incoming atom and       |"
  write(*,*) "setting up simulation box...                       |"
  write(*,*) "                                                   |"
! The box should be "padding" Ansgt taller than the highest z coordinate
! and "padding" Angst shorter than the shortest z coordinate
  lz(2) = -1.d10
  lz(1) = 1.d10
  do i = 1, n_atoms
    if( .not. remove_atom(i) )then
      if( positions(i, 3) > lz(2) )then
        lz(2) = positions(i, 3)
      end if
      if( positions(i, 3) < lz(1) )then
        lz(1) = positions(i, 3)
      end if
    end if
  end do
  lz(2) = lz(2) + padding
  lz(1) = lz(1) - padding


  if( add_ion )then
!   We calculate the incident ion's properties
    vel(1:2) = 0.d0
    vel(3) = -dsqrt(2.d0*ion_energy*1.602/ion_mass/1.66)*1.d2
!   Seed random number generator
    call random_seed(size = n)
    allocate( seed(n) )
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    call random_number(rand)
    pos(1) = lx(1) + rand*(lx(2) - lx(1))
    call random_number(rand)
    pos(2) = ly(1) + rand*(ly(2) - ly(1))
!   Do not waste time. Place the incident ion min_dist Angst away from any neighboring atom
    pos(3) = 0.d0
    do i = 1, n_atoms
      call get_distance( (/ pos(1), pos(2), 0.d0 /), (/ positions(i,1), positions(i,2), 0.d0 /), &
                          lx, ly, lz, dxy )
      if( dxy < min_dist )then
        dz = dsqrt( min_dist**2 - dxy**2 )
        if( positions(i, 3) + dz > pos(3) )then
          pos(3) = positions(i, 3) + dz
        end if
      end if
    end do
  end if

! Print the new structure
  open(unit=10, file=input, status="unknown")
  if( add_ion )then
    n = n_atoms - n_remove + 1
  else
    n = n_atoms - n_remove
  end if
  write(10,*)
  write(10,*) n, "atoms"
  write(10,*)
  write(10,*) "1 atom types"
  write(10,*)
  write(10,*) lx(1:2), "xlo xhi"
  write(10,*) ly(1:2), "ylo yhi"
  write(10,*) lz(1:2), "zlo zhi"
  write(10,*)
  write(10,*) "Masses"
  write(10,*)
  write(10,*) "1", ion_mass
  write(10,*)
  write(10,*) "Atoms"
  write(10,*)
  j = 0
  do i = 1, n_atoms
    if( .not. remove_atom(i) )then
      j = j + 1
      write(10,*) j, "1", positions(i, 1:3)
    end if
  end do
  if( add_ion )then
    write(10,*) n, "1", pos(1:3)
  end if
  write(10,*)
  write(10,*) "Velocities"
  write(10,*)
  j = 0
  do i = 1, n_atoms
    if( .not. remove_atom(i) )then
      j = j + 1
      write(10,*) j, velocities(i, 1:3)
    end if
  end do
  if( add_ion )then
    write(10,*) n, vel(1:3)
  end if
  close(10)


  write(*,*) "***************************************************|"
  write(*,*) "                                                   |"
  write(*,*) "Program successfully executed.                     |"
  write(*,*) "                                                   |"
  write(*,*) "___________________________________________________/"


end program
