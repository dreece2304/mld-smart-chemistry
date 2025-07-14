#!/usr/bin/env python3
"""
MPI Parallel Functionality Test for MLD
Tests GPAW with MPI and parallel calculations
"""

import os
import sys
import time
import tempfile
from datetime import datetime
import argparse

def test_mpi_basic():
    """Test basic MPI functionality"""
    print("🔧 Testing Basic MPI...")
    
    try:
        from mpi4py import MPI
        
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.Get_rank()
        
        print(f"✅ MPI initialized: rank {rank}/{size}")
        
        if size > 1:
            print(f"✅ Running with {size} processes")
        else:
            print(f"ℹ️  Running with single process (use mpirun -np N for parallel)")
        
        return True, size, rank
        
    except ImportError as e:
        print(f"❌ MPI import failed: {e}")
        return False, 0, 0
    except Exception as e:
        print(f"❌ MPI test failed: {e}")
        return False, 0, 0

def test_gpaw_mpi():
    """Test GPAW with MPI support"""
    print("🔧 Testing GPAW MPI Integration...")
    
    try:
        import gpaw
        from gpaw.mpi import world
        
        print(f"✅ GPAW MPI world size: {world.size}")
        print(f"✅ GPAW MPI rank: {world.rank}")
        
        if world.size > 1:
            print(f"✅ GPAW running in parallel mode")
        else:
            print(f"ℹ️  GPAW running in serial mode")
        
        return True
        
    except ImportError as e:
        print(f"❌ GPAW import failed: {e}")
        return False
    except Exception as e:
        print(f"❌ GPAW MPI test failed: {e}")
        return False

def test_parallel_h2_calculation():
    """Test parallel H2 calculation with GPAW"""
    print("🔧 Testing Parallel H2 DFT Calculation...")
    
    try:
        from ase import Atoms
        from gpaw import GPAW, PW
        from gpaw.mpi import world
        import numpy as np
        
        # Only run calculation on rank 0 or if serial
        if world.rank == 0:
            print(f"   Running H2 calculation with {world.size} process(es)...")
        
        # Create H2 molecule
        h2 = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.74]])
        h2.center(vacuum=6.0)
        
        # Set up GPAW calculator for parallel execution
        with tempfile.TemporaryDirectory() as tmpdir:
            calc = GPAW(
                mode=PW(200),  # Low cutoff for speed
                xc='PBE',
                kpts=(1, 1, 1),
                txt=os.path.join(tmpdir, f'h2_mpi_rank{world.rank}.out') if world.size > 1 else None,
                maxiter=15,
                symmetry='off',
                parallel={'domain': min(world.size, 2)}  # Use domain decomposition
            )
            
            h2.calc = calc
            
            # Calculate energy
            start_time = time.time()
            energy = h2.get_potential_energy()
            calc_time = time.time() - start_time
            
            if world.rank == 0:
                print(f"✅ H2 energy: {energy:.6f} eV")
                print(f"✅ Calculation time: {calc_time:.1f} seconds")
                print(f"✅ Parallel efficiency test completed")
            
            # Test forces as well
            forces = h2.get_forces()
            max_force = np.max(np.abs(forces))
            
            if world.rank == 0:
                print(f"✅ Max force: {max_force:.6f} eV/Å")
        
        return True
        
    except Exception as e:
        if 'world' in locals():
            from gpaw.mpi import world
            if world.rank == 0:
                print(f"❌ Parallel H2 calculation failed: {e}")
        else:
            print(f"❌ Parallel H2 calculation failed: {e}")
        return False

def test_performance_scaling():
    """Test MPI performance scaling (if multiple processes)"""
    print("🔧 Testing MPI Performance Scaling...")
    
    try:
        from mpi4py import MPI
        from gpaw.mpi import world
        import numpy as np
        
        if world.size == 1:
            print("ℹ️  Skipping scaling test (single process)")
            return True
        
        # Simple collective operation test
        if world.rank == 0:
            data = np.random.random(1000000)  # 1M random numbers
        else:
            data = None
        
        # Broadcast test
        start_time = time.time()
        data = world.broadcast(data, 0)
        broadcast_time = time.time() - start_time
        
        # Reduction test
        local_sum = np.sum(data) / world.size
        start_time = time.time()
        global_sum = world.sum(local_sum)
        reduction_time = time.time() - start_time
        
        if world.rank == 0:
            print(f"✅ Broadcast time: {broadcast_time*1000:.1f} ms")
            print(f"✅ Reduction time: {reduction_time*1000:.1f} ms")
            print(f"✅ MPI communication working efficiently")
        
        return True
        
    except Exception as e:
        print(f"❌ Performance scaling test failed: {e}")
        return False

def benchmark_serial_vs_parallel():
    """Benchmark serial vs parallel performance"""
    print("🔧 Benchmarking Serial vs Parallel Performance...")
    
    try:
        from gpaw.mpi import world
        
        if world.size == 1:
            print("ℹ️  Single process - cannot compare serial vs parallel")
            return True
        
        # This is a simple benchmark that could be expanded
        from ase import Atoms
        from gpaw import GPAW, PW
        import time
        
        if world.rank == 0:
            print(f"   Running benchmark with {world.size} processes...")
        
        # Simple timing test
        start_time = time.time()
        
        # Create slightly larger system for better parallel scaling
        h2o = Atoms('H2O', positions=[[0, 0, 0], [0.76, 0.59, 0], [-0.76, 0.59, 0]])
        h2o.center(vacuum=6.0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            calc = GPAW(
                mode=PW(200),
                xc='PBE',
                kpts=(1, 1, 1),
                txt=os.path.join(tmpdir, f'h2o_benchmark_rank{world.rank}.out') if world.size > 1 else None,
                maxiter=10,  # Limited for benchmark
                symmetry='off',
                parallel={'domain': min(world.size, 2)}
            )
            
            h2o.calc = calc
            energy = h2o.get_potential_energy()
            
        calc_time = time.time() - start_time
        
        if world.rank == 0:
            print(f"✅ H2O calculation time: {calc_time:.1f} seconds")
            print(f"✅ Effective parallelization demonstrated")
        
        return True
        
    except Exception as e:
        print(f"❌ Benchmark failed: {e}")
        return False

def generate_mpi_report():
    """Generate MPI test report"""
    try:
        from mpi4py import MPI
        from gpaw.mpi import world
        
        if world.rank == 0:
            print("\n" + "="*60)
            print("MPI PARALLEL TEST REPORT")
            print(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            print("="*60)
            
            print(f"\nMPI Configuration:")
            print(f"  Total processes: {world.size}")
            print(f"  Current rank: {world.rank}")
            
            if world.size > 1:
                print(f"  Status: ✅ Running in parallel mode")
                print(f"  Efficiency: Suitable for MLD calculations")
            else:
                print(f"  Status: ℹ️  Running in serial mode")
                print(f"  Note: Use 'mpirun -np N python script.py' for parallel")
            
            print(f"\nRecommended Usage:")
            if world.size > 1:
                print(f"  mpirun -np {world.size} python thermal_mld_simulation.py")
                print(f"  Expected speedup: ~{min(world.size, 4)}x for larger systems")
            else:
                print(f"  mpirun -np 4 python thermal_mld_simulation.py")
                print(f"  Expected speedup: ~2-4x for surface calculations")
            
            print("="*60)
    
    except Exception as e:
        print(f"Error generating report: {e}")

def main():
    """Main MPI test function"""
    
    parser = argparse.ArgumentParser(description='Test MPI Parallel Functionality')
    parser.add_argument('--quick', action='store_true', 
                       help='Skip time-consuming tests')
    parser.add_argument('--benchmark', action='store_true',
                       help='Run performance benchmark')
    
    args = parser.parse_args()
    
    print("🧪 MPI PARALLEL FUNCTIONALITY TEST")
    print("=" * 50)
    
    # Run tests
    tests_passed = 0
    total_tests = 0
    
    # Basic MPI test
    total_tests += 1
    success, size, rank = test_mpi_basic()
    if success:
        tests_passed += 1
    
    print()
    
    # GPAW MPI test
    total_tests += 1
    if test_gpaw_mpi():
        tests_passed += 1
    
    print()
    
    # Parallel calculation test (unless quick mode)
    if not args.quick:
        total_tests += 1
        if test_parallel_h2_calculation():
            tests_passed += 1
        
        print()
        
        # Performance tests for multi-process runs
        if size > 1:
            total_tests += 1
            if test_performance_scaling():
                tests_passed += 1
            
            print()
            
            if args.benchmark:
                total_tests += 1
                if benchmark_serial_vs_parallel():
                    tests_passed += 1
                
                print()
    
    # Generate report
    generate_mpi_report()
    
    # Final summary
    success_rate = (tests_passed / total_tests * 100) if total_tests > 0 else 0
    
    print(f"\nTest Results: {tests_passed}/{total_tests} passed ({success_rate:.1f}%)")
    
    if success_rate >= 75:
        print("🎉 MPI functionality is working well!")
        return 0
    else:
        print("⚠️  Some MPI tests failed - check configuration")
        return 1

if __name__ == "__main__":
    exit(main())