#include "X86.h"
#include "X86InstrInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include "llvm/CodeGen/MachineRegisterInfo.h"
#include "llvm/CodeGen/MachineFrameInfo.h"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <float.h>

using namespace llvm;

#define X86_REGISTER_ALLOCATOR_PASS_DESC "A basic register allocator for X86"
#define X86_REGISTER_ALLOCATOR_PASS_NAME "x86-reg-alloc"

namespace {

enum EXIT_STATUS_T {DONE, FAILED, SPILLED};

class LiveRange {
public:
    // std::vector<MachineBasicBlock*> BBs;
    // list of instructions where this equivalence class is live-in
    std::set<MachineInstr*> liveInstrs;

    // regs in one  equivalence class
    std::set<unsigned> virtRegs;
    
    // the physical register allocated to this LiveRange
    unsigned physicalReg;
    unsigned degree;
    bool isAllocated;
    bool usesPhysicalReg;

    LiveRange() {usesPhysicalReg = false;isAllocated = false; degree = 0;}
    LiveRange(unsigned virtReg) : LiveRange() {
        virtRegs.insert(virtReg);
    }
    ~LiveRange() {}


    bool contains (unsigned virtReg) {
        return virtRegs.count(virtReg);
    }

    bool contains (std::set<unsigned>& _virtRegs) {
        auto __first1 = virtRegs.begin();
        auto __first2 = _virtRegs.begin();

        while (__first1 != virtRegs.end() && __first2 != _virtRegs.end())
            if (*__first1 < *__first2)
              ++__first1;
            else if (*__first2 < *__first1)
              ++__first2;
            else
              return true;

        return false;
    }

    bool contains (MachineInstr* MI) {
        return liveInstrs.count(MI);
    }    
    
    double calculateAndGetCostRatio(MachineRegisterInfo& regInfo) {
    // formula used = #defs+uses / degree 
        double useDefCount = 0;
        for (unsigned virtReg : virtRegs) {
            for (auto& useDef : regInfo.reg_operands(virtReg)) {
                (void)useDef;
                ++useDefCount;
            }
        }
        return (useDefCount/degree);
    }

    void MergeWith (LiveRange* LR) {
        liveInstrs.insert(LR->liveInstrs.begin(),LR->liveInstrs.end());
        virtRegs.insert(LR->virtRegs.begin(),LR->virtRegs.end());
    }

    bool interferesWith (LiveRange* LR) {
        // doing intersection is a waste of time
        // replaced with a search that stops at first matching element
        // sets are ordered

        std::vector<MachineInstr*> intersectionSet;
        // std::set_intersection (liveInstrs.begin(), liveInstrs.end(), 
            // LR->liveInstrs.begin(), LR->liveInstrs.end(), intersectionSet.begin());

        auto __first1 = liveInstrs.begin();
        auto __first2 = LR->liveInstrs.begin();

      while (__first1 != liveInstrs.end() && __first2 != LR->liveInstrs.end())
        if (*__first1 < *__first2)
          ++__first1;
        else if (*__first2 < *__first1)
          ++__first2;
        else
          return true;

      return false;
        // return (intersectionSet.size() > 0);
    }

    void addLiveInstr (MachineInstr* MI) {
        liveInstrs.insert(MI);
    }

    void addLiveInstrs (std::set<MachineInstr*> _liveInstrs) {
        liveInstrs.insert(_liveInstrs.begin(), _liveInstrs.end());
    }
};


class X86RegisterAllocator : public MachineFunctionPass {
public:
    static char ID;

    // might want to change to a vector
    std::set<LiveRange*> interferenceGraph;
    std::map<unsigned, LiveRange*> interferenceGraphRetrieval;
    std::map<MachineInstr*, std::set<unsigned>> livenessInformation;
    std::set<unsigned> spilledRegisters;

    // list of phi nodes that we need to remove after reg allocation
    std::set<MachineInstr*> phiNodes;

    X86RegisterAllocator() : MachineFunctionPass(ID) {
        initializeX86RegisterAllocatorPass(*PassRegistry::getPassRegistry());
    }

    StringRef getPassName() const override { 
        return X86_REGISTER_ALLOCATOR_PASS_NAME; 
    }

    bool runOnMachineFunction(MachineFunction &MF) override;

    EXIT_STATUS_T color(MachineFunction &MF, MachineRegisterInfo& regInfo);

    void calcGlobalLivenessInfo(MachineFunction &MF, MachineRegisterInfo& regInfo);

    void calcLivenessForDef(MachineInstr* def, MachineOperand& defOp, MachineRegisterInfo& regInfo);

    std::set<MachineInstr*> getPhiNodes() {
        return phiNodes;
    }
    
    void addPhiNode (MachineInstr* MI) {
        phiNodes.insert(MI);
    }

    LiveRange* makeNewLiveRangeFor(unsigned _physicalReg) {
        LiveRange* LR = new LiveRange();
        interferenceGraph.insert(LR);
        LR->physicalReg = _physicalReg;
        LR->usesPhysicalReg = true;
        return LR;
    }

    LiveRange* getLiveRangeFor(unsigned virtReg) {
        if (interferenceGraphRetrieval.find(virtReg) != interferenceGraphRetrieval.end()) {
            return interferenceGraphRetrieval[virtReg];
        } else {
            // live range does not exist so create new one
            LiveRange* LR = new LiveRange(virtReg);
            interferenceGraph.insert(LR);
            interferenceGraphRetrieval[virtReg] = LR;

            return LR;
        }
    }

    std::vector<LiveRange*> getLiveRangesFor(MachineInstr* MI) {
        std::vector<LiveRange*> ranges;
        for (auto& LR :interferenceGraph) {
            if (LR->contains(MI)) {
                ranges.push_back(LR);
            }
        }

        return ranges;
    }

    std::vector<LiveRange*> getInterferingLiveRangesFor(LiveRange* LR) {
        std::vector<LiveRange*> interferingRanges;
        for (auto& otherLR :interferenceGraph) {
            if (otherLR == LR) continue;
            if (otherLR->interferesWith(LR)) {
                interferingRanges.push_back(otherLR);
            }
        }

        return interferingRanges;
    }

    MachineBasicBlock::iterator getIterForMI(MachineInstr& MI) {
        auto MBB = MI.getParent();
        MachineBasicBlock::iterator i = MBB->begin();
        for (; i != MBB->end(); ++i) {
            if (MI.isIdenticalTo(*i)) break;  // will always break
        }
        return i;
    }

};

bool X86RegisterAllocator::runOnMachineFunction(MachineFunction &MF) {
    outs() << "\n" << "=====================\n" << MF.getName() << "()\n=====================" <<"\n";
    for (auto &MBB : MF) {
        outs() << "Contents of MachineBasicBlock:\n";
        outs() << MBB << "\n";
        outs().flush();
    }
    // cleanup from previous function run...if any
    spilledRegisters.clear();
BEGINNING:
    for (LiveRange* LR : interferenceGraph) {
        delete LR;
    }
    interferenceGraph.clear();
    interferenceGraphRetrieval.clear();
    livenessInformation.clear();
    phiNodes.clear();

    MachineRegisterInfo& regInfo = MF.getRegInfo();
    calcGlobalLivenessInfo(MF, regInfo);

    if (interferenceGraph.size() > 0 && color(MF, regInfo) == EXIT_STATUS_T::SPILLED) {
        // return true;
        goto BEGINNING;
    }

    for (MachineInstr* MI : getPhiNodes()) {
        MI->eraseFromParent();
    }
    MF.getProperties().set(MachineFunctionProperties::Property::NoPHIs);
    MF.getProperties().set(MachineFunctionProperties::Property::NoVRegs);

    outs() << "\n\n";
    for (auto &MBB : MF) {
        outs() << "Contents of MachineBasicBlock:\n";
        outs() << MBB << "\n";
        outs().flush();
    }
    return true;
}


EXIT_STATUS_T X86RegisterAllocator::color(MachineFunction &MF, MachineRegisterInfo& regInfo) {

    const TargetRegisterInfo* tRegInfo = MF.getSubtarget().getRegisterInfo();
    const TargetInstrInfo* tInstrInfo = MF.getSubtarget().getInstrInfo();
    
    // get any virtual reg num; just for the purpose of getting the allocatable reg list size
    // physical only liveranges wont have any virtRegs though
    // so we have to iterate
    unsigned temp_virtReg;
    bool found_Temp_virtReg = false;
    for (LiveRange* LR : interferenceGraph) {
        if (!LR->usesPhysicalReg) {
        // if (LR->virtRegs.size() > 0) {
            found_Temp_virtReg = true;
            temp_virtReg = *(LR->virtRegs.begin());
            break;
        }
    }

    if (!found_Temp_virtReg) {
        // means that there are no virtual ranges! 
        // nothing to do!!
        return EXIT_STATUS_T::DONE;
    }
    // outs() << "temp_virtReg: " << TargetRegisterInfo::virtReg2Index(temp_virtReg) << "\n";
    BitVector temp_regBitVec = tRegInfo->getAllocatableSet(MF, regInfo.getRegClass(temp_virtReg));
    unsigned numPhysRegs = temp_regBitVec.count();

    // misc code for checking number of registers in different register classes
    // all of them (64-bit, 32-bit, 8-bit) are returning 15 for my sample code

    // for (std::set<LiveRange*>::iterator i = interferenceGraph.begin(); i != interferenceGraph.end();++i) {
    //         for (unsigned virtReg : (*i)->virtRegs) {
    //             const TargetRegisterClass* tRegClass = regInfo.getRegClass(virtReg);
    //             BitVector temp_regBitVec = MF.getSubtarget().getRegisterInfo()->getAllocatableSet(MF, tRegClass);
    //              numPhysRegs = temp_regBitVec.count();
    // outs() << "Total of N = " << numPhysRegs << " Registers in architecture\n";
    //     for (auto index : temp_regBitVec.set_bits()) {
    //     outs() << "Reg: " << MF.getSubtarget().getRegisterInfo()->getRegAsmName(index) << "\n";
    // }
    //         }
    //     }

    outs() << "*** Total of N = " << numPhysRegs << " Allocatable Registers in architecture\n";
    for (auto index : temp_regBitVec.set_bits()) {
        outs() << "\tAllocatable Reg: " << tRegInfo->getRegAsmName(index) << "\n";
    }
    //// now print regs that are not Allocatable
    // temp_regBitVec = temp_regBitVec.flip();

    // for (auto index : temp_regBitVec.set_bits()) {
    //     outs() << "\tNOT-Allocatable Reg: " << MF.getSubtarget().getRegisterInfo()->getRegAsmName(index) << "\n";
    // }


    std::stack<LiveRange*> colorStack;

    // too time consuming
    // todo: convert to a 2d matrix and update when necessary
    outs() << "*** Building the interference graph\n";
    for (std::set<LiveRange*>::iterator i = interferenceGraph.begin(); i != interferenceGraph.end(); ++i) {
        for (std::set<LiveRange*>::iterator j = i; j != interferenceGraph.end(); ++j) {
            if (i==j) continue;

            if ((*i)->interferesWith((*j))) {
                ((*i)->degree)++;
                ((*j)->degree)++;
            }
        }
    }

    outs() << "*** Pushing nodes to the stack for the colouring. ";
    outs() << "Total LiveRanges: " << interferenceGraph.size() << "\n";
    
    while (interferenceGraph.size() > 0) {
        bool removedSomething = false;

        for (std::set<LiveRange*>::iterator i = interferenceGraph.begin(); i != interferenceGraph.end();) {
            if ((*i)->degree < numPhysRegs) {
                removedSomething = true;
                colorStack.push((*i));

                // technically we have to remove the node and its edges but
                // we do not have an explicit/traditional `graph` with edges
                // so we just decrement the degree for the interfering ranges
                i = interferenceGraph.erase(i);
                for (std::set<LiveRange*>::iterator j = interferenceGraph.begin(); j != interferenceGraph.end();++j) {
                    if ((*j)->interferesWith(colorStack.top())) {
                        ((*j)->degree)--;
                    }
                }
            } else {
                ++i;
            }
        }

        if (!removedSomething) {
            outs() << "\t*** Could not find any node with degree < " << numPhysRegs << "...trying to spill\n";  
        } else if (interferenceGraph.size() > 0) {
            outs() << "\t*** Removed some nodes...now trying to push nodes from the reduced interference graph. ";
            outs() << "Total LiveRanges: " << interferenceGraph.size() << "\n";
        } else {
            outs() << "\t*** Removed some nodes...All nodes pushed to the stack!\n";
        }
 
        if (!removedSomething) {
            // calculate cost:benefit ratio
            // pick node with lowest number
            // formula used = #defs+uses / degree 
            // todo also look at nest depth 
            LiveRange* spillRange = NULL;
            double costRatio = DBL_MAX;
            double virtualRangesExist = false;

            for (LiveRange* LR : interferenceGraph) {
                if (LR->usesPhysicalReg == false) {
                    virtualRangesExist = true;
                    break;
                }
            }

            outs() << "\t\t*** Calculating the cost/benefit ratio\n";
            for (std::set<LiveRange*>::iterator i = interferenceGraph.begin(); i != interferenceGraph.end();++i) {
                // dont want to choose already spilled ranges 
                // eventhough they might have the lowest cost
                if ((*i)->contains(spilledRegisters)) {
                    continue;
                }
                // dont spill physical ranges
                if (virtualRangesExist && (*i)->usesPhysicalReg) {
                    continue;
                }
                double currCost = (*i)->calculateAndGetCostRatio(regInfo);
                if (currCost < costRatio) {
                    spillRange = *i;
                    costRatio = currCost;
                }
            }

            if (spillRange == NULL)
                report_fatal_error("FAILED!! Could not find any node to spill");
            // assert(spillRange != NULL && "Could not find any node to spill");
            // if (spillRange == NULL) {
            //     return EXIT_STATUS_T::FAILED;
            // }

            outs() << "\t\t*** Spilling live range class for reg #'s: ";
            for (unsigned virtReg : spillRange->virtRegs) {
                outs() << TargetRegisterInfo::virtReg2Index(virtReg) << " ";
            }
            outs() << ". Degree = " << spillRange->degree << "\n";

            // todo also check if we are picking something that is already on the stack
            unsigned _virtReg = *(spillRange->virtRegs.begin());
            const TargetRegisterClass* tRegClass = 
            spillRange->usesPhysicalReg? tRegInfo->getMinimalPhysRegClass(spillRange->physicalReg) : regInfo.getRegClass(_virtReg);

            int stackSlot = MF.getFrameInfo().CreateSpillStackObject(
                                tRegInfo->getSpillSize(*tRegClass), 
                                tRegInfo->getSpillAlignment(*tRegClass));

            // for all virtual registers in the liverange equivalence class
            // for all defining instructions: save the definition into the same stackSlot
            // for all uses: add load instructions from the stack
            // ignore phi nodes; delete them later
            for (unsigned virtReg : spillRange->virtRegs) {
                for (MachineInstr& MI : regInfo.use_instructions(virtReg)) {
                    // todo check if deleting uses in a range loop causes any problem
                        // ignore if phi
                    if (!MI.isPHI()) {
                        unsigned newVirtReg = regInfo.createVirtualRegister(tRegClass);
                        spilledRegisters.insert(newVirtReg);
                        // load instruction is added BEFORE the specified machine instruction
                        // we need to add a load before the use
                        tInstrInfo->loadRegFromStackSlot (
                            /*MachineBasicBlock&*/ *(MI.getParent()),
                            /*MachineBasicBlock::iterator*/ getIterForMI(MI),
                            /*DestReg*/ newVirtReg,
                            /*FrameIndex*/ stackSlot,
                            /*const TargetRegisterClass**/     tRegClass,
                            /*const TargetRegisterInfo**/      tRegInfo);

                        // replace this use's operand(s) by this load
                        for (MachineOperand& use : MI.uses()) {
                            if (!use.isReg()) continue;
                            if (!TargetRegisterInfo::isVirtualRegister(use.getReg())) continue;
                            if (use.getReg() == virtReg) {
                                use.setReg(newVirtReg);
                            }
                        }            
                    } else {
                        phiNodes.erase(&MI);

                        MI.eraseFromParent();
                    }
                }  // end loop for spilling uses

                // defs come later because deleting PHI nodes prematurely might
                // mess with its uses...also we kill the register

                for (MachineInstr& MI : regInfo.def_instructions(virtReg)) {
                    if (MI.isPHI()) {
                        // no longer need a phi node since all uses and defs share a stack slot
                        phiNodes.erase(&MI);
                        
                        // todo will this break the for loop?
                        MI.eraseFromParent();
                    } else {
                        // store instruction is added BEFORE the specified machine instruction
                        // we need to add a store AFTER each def
                        MachineBasicBlock::iterator it = getIterForMI(MI);
                        ++it;

                        spilledRegisters.insert(virtReg);
                        tInstrInfo->storeRegToStackSlot (
                            /*MachineBasicBlock&*/ *(MI.getParent()),
                            /*MachineBasicBlock::iterator*/ it,
                            /*SrcReg*/ virtReg,
                            /*isKill*/ false,
                            /*FrameIndex*/ stackSlot,
                            /*const TargetRegisterClass**/     tRegClass,
                            /*const TargetRegisterInfo**/      tRegInfo);
                    }
                }  // end loop for spilling uses

            }  // end loop for spilling all virtRegs in the equivalence class

            // now that we have spilled one live range class we need to restart the algorithm
            return EXIT_STATUS_T::SPILLED;
        }  // end if (!removedSomething) in the stack pushing phase
    }  // end while (interferenceGraph.size() > 0)  --> the stack pushing phase

    // now start popping from the stack and assign actual physical registers
    // first remove all the physical ranges since they already have an assigned register
    // we cant iterate on a stack so rebuild the structures

    std::stack<LiveRange*> colorStack2;
    while (!colorStack.empty()) {
        LiveRange* LR = colorStack.top();
        colorStack.pop();

        if (LR->usesPhysicalReg) {
            LR->isAllocated = true;
            // rebuild the interferenceGraph
            interferenceGraph.insert(LR);
        } else {
            colorStack2.push(LR);
        }
    }

    while (!colorStack2.empty()) {
        colorStack.push(colorStack2.top());
        colorStack2.pop();
    }

    while (!colorStack.empty()) {
        LiveRange* LR = colorStack.top();
        colorStack.pop();
        unsigned currVirtReg = *(LR->virtRegs.begin());
        // auto tRegClass = regInfo.getRegClass(currVirtReg);
        // if (tRegClass->isASubClass()) {
        // iterator sorted by order so we get the biggest class
            // outs() << tRegInfo->getRegClassName (tRegClass) << "\n";
            // outs() << tRegInfo->getRegClassName (tRegInfo->getLargestLegalSuperClass (tRegClass, MF)) << "\n";

            // for (auto i = tRegClass->getSuperClasses(); *i != NULL; ++i) {
                // outs() << tRegInfo->getRegClassName (*i) << "\n";
            // tRegClass = *i;
            // }
            // tRegClass = tRegInfo->getLargestLegalSuperClass (tRegClass, MF);
            // tRegClass = *(tRegClass->getSuperClasses());
        // }

        // make sure all instruction contraints are taken into account before getting allocatable registers
        // first constrain the first virtual register
        for (unsigned virtReg : LR->virtRegs) {
            regInfo.constrainRegAttrs(currVirtReg, virtReg);
        }
        // then constrain all other virtual registers 
        // with the first register's contraints
        for (unsigned virtReg : LR->virtRegs) {
            regInfo.constrainRegAttrs(virtReg, currVirtReg);
        }

        BitVector temp_regBitVec = tRegInfo->getAllocatableSet(MF, 
            LR->usesPhysicalReg? 
                tRegInfo->getMinimalPhysRegClass(LR->physicalReg) : regInfo.getRegClass(currVirtReg)
            );
        
        // save the allocatable physical registers in a set for easy comparison
        std::set<unsigned> allocatableRegisters;
        std::set<unsigned> unAllocatableRegisters;
        for (auto index : temp_regBitVec.set_bits()) {
            allocatableRegisters.insert(index);
        }


        for (LiveRange* otherLR : interferenceGraph) {
            if (LR->interferesWith(otherLR)) {
                // outs() << TargetRegisterInfo::virtReg2Index(*(LR->virtRegs.begin())) << " interferesWith " << TargetRegisterInfo::virtReg2Index(*(otherLR->virtRegs.begin())) << "\n";

                // find which register it overlaps with and mark it as unAllocatable
                if (allocatableRegisters.count(otherLR->physicalReg)) {
                    unAllocatableRegisters.insert(otherLR->physicalReg);
                } else {
                    for (unsigned r : allocatableRegisters) {
                        if (tRegInfo->regsOverlap(r, otherLR->physicalReg)) {
                            unAllocatableRegisters.insert(r);
                            // break;
                        }
                    }
                }
            }
        }

        // give all the virtual registers in the live range
        // any remaining color from the allocatableRegisters
        for (unsigned r : unAllocatableRegisters) {
            allocatableRegisters.erase(r);
        }

        if (LR->usesPhysicalReg) {
            if (!allocatableRegisters.count(LR->physicalReg)) {
                outs() << "Uses: " << tRegInfo->getRegAsmName(LR->physicalReg) << "\n";
                for (unsigned a : allocatableRegisters)
                    outs() << "allocatableRegisters Left: " << tRegInfo->getRegAsmName(a) << "\n";
                report_fatal_error("Could not color.");
            }
        } else {
            LR->physicalReg = *(allocatableRegisters.begin());
            outs() << "*** Allocated: " << tRegInfo->getRegAsmName(LR->physicalReg) << "\n";
        }
        LR->isAllocated = true;

        for (unsigned virtReg : LR->virtRegs) {
            regInfo.replaceRegWith(virtReg, LR->physicalReg);   
        }

        // rebuild the interferenceGraph
        interferenceGraph.insert(LR);
    }

    return EXIT_STATUS_T::DONE;
}

void X86RegisterAllocator::calcGlobalLivenessInfo(MachineFunction &MF, MachineRegisterInfo& regInfo) {
    
    // liveness analysis for every def/virtual SSA Reg 
    for (auto &MBB : MF) {
        // entry block has physical registers which are live in due to x86 calling convention
        // they are copied to a virtual reg asap
        // assume killed afterwards
        if (MBB.pred_empty()) {  // entry block
            for (/*RegisterMaskPair*/auto regPair : MBB.liveins()) {
                auto regNum = regPair.PhysReg;
                std::set<MachineInstr*> liveInstrs;
                LiveRange* LR = makeNewLiveRangeFor(regNum);

                for (auto& MI : MBB) {
                    liveInstrs.insert(&MI);
                        
                    bool foundUse = false;
                    for (MachineOperand& u : MI.uses()) {
                        if (!u.isReg()) continue;
                        if (!TargetRegisterInfo::isPhysicalRegister(u.getReg())) continue;
                        if (u.getReg() == regNum) {
                            foundUse = true;
                            break;
                        }
                    }   
                    if (foundUse) break;                 
                }
                // update the live range
                LR->addLiveInstrs(liveInstrs);

                // update the MI-to-virtRegs liveness map
                for (MachineInstr* MI : liveInstrs) {
                    livenessInformation[MI].insert(regNum);
                }
            }
        }

        for (auto& MI : MBB) {
            // only gives explicit defs
            // e.g. move into args before function call
            for (MachineOperand& def : MI.defs()) {
                if (!def.isReg()) continue;
                if (TargetRegisterInfo::isPhysicalRegister(def.getReg())) {
                    // for physical regs...case where there are defs
                    // live range is only from the def until the next use within the BB
                    std::set<MachineInstr*> liveInstrs;
                    LiveRange* LR = makeNewLiveRangeFor(def.getReg());

                    auto i = getIterForMI(MI);
                    bool foundUse = false;
                    // some defs are dead e.g. MUL instruction
                    if (def.isDead()) {
                        liveInstrs.insert(&(*i)); 
                        foundUse = true;
                    }
                        
                    for (++i; i!=MBB.end(); ++i) {
                        if (foundUse) break;
                        liveInstrs.insert(&(*i));
                        // also includes implicit uses
                        for (MachineOperand& u : i->uses()) {
                            if (!u.isReg()) continue;
                            if (!TargetRegisterInfo::isPhysicalRegister(u.getReg())) continue;
                            if (u.getReg() == def.getReg()) {
                                // LR->virtRegs.insert(i->operands_begin()->getReg());
                                foundUse = true;
                                break;
                            }
                        }
                    }
                    // update the live range
                    LR->addLiveInstrs(liveInstrs);

                    // update the MI-to-virtRegs liveness map
                    for (MachineInstr* MI : liveInstrs) {
                        livenessInformation[MI].insert(def.getReg());
                    }
                }
                if (!TargetRegisterInfo::isVirtualRegister(def.getReg())) continue;
                    
                // perform Liveness analysis and 
                // update the LiveRange for the virtReg, and
                // update the MI to virtRegs liveness map
                // unsigned regNum = def.getReg();
                // outs() << "*** calcLivenessForDef reg #" << TargetRegisterInfo::virtReg2Index(regNum) << "\n";

                calcLivenessForDef(&MI, def, regInfo);
            }

            // also for physical registers
            // the case above does not look at implicit definitions
            // physical reg %rax, used for returns is an implicit definition
            for (MachineOperand& MO : MI.implicit_operands()) {
                if (MO.isReg() && MO.isDef() && TargetRegisterInfo::isPhysicalRegister(MO.getReg())) {
                    // not are implicit defs are allocatable
                    // e.g. %rsp etc
                    if (regInfo.isAllocatable(MO.getReg())) {
                        // live range is the call site, 
                        // and, if ret val is used, all instructions until that one 
                        std::set<MachineInstr*> liveInstrs;
                        LiveRange* LR = makeNewLiveRangeFor(MO.getReg());

                        auto i = getIterForMI(MI);
                        liveInstrs.insert(&(*i));

                        bool foundUse = false;
                        if (MO.isDead()) {
                            foundUse = true;
                        }
                            
                        for (++i; i!=MBB.end(); ++i) {
                            if (foundUse) break;
                            liveInstrs.insert(&(*i));
                            // also includes implicit uses
                            for (MachineOperand& u : i->uses()) {
                                if (!u.isReg()) continue;
                                if (!TargetRegisterInfo::isPhysicalRegister(u.getReg())) continue;
                                if (u.getReg() == MO.getReg()) {
                                    // LR->virtRegs.insert(i->operands_begin()->getReg());
                                    foundUse = true;
                                    break;
                                }
                            }
                        }

                        // update the live range
                        LR->addLiveInstrs(liveInstrs);
            // outs() << "PHYS:: retval implicit defs num instructions: " <<  LR->liveInstrs.size() << "\n\n";

                        // update the MI-to-virtRegs liveness map
                        for (MachineInstr* MI : liveInstrs) {
                            livenessInformation[MI].insert(MO.getReg());
                        }                
                    }
                }
            }
                
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////
    // for (MachineInstr* phi : getPhiNodes()) {
        // outs() << "PHI: "<< *phi << "\n";
        // outs() << "getNumOperands: "<< phi->getNumOperands() << "\n";
        // for () {
        //     outs() << "operand "<<i<<": "<< phi->getOperand(i) << "\n";
        // }

        /*
            outs() << "PHI: "<< *user << "\n";
            outs() << "index of reg is: " << user->findRegisterUseOperandIdx(regNum) << "\n"; 
            outs() << "BB is of reg is: " <<  << "\n"; 
            // PHI.getNumOperands()
// return PHI.getOperand(Index * 2 + 2).getMBB();

        */
    //     for (MachineOperand& def : phi->defs()) {
    //         if (!def.isReg()) continue;
    //         if (!TargetRegisterInfo::isVirtualRegister(def.getReg())) continue;
    //         outs() << "\tPHI: defines: "<< def << "\n";
    //     }

    //     for (auto& use : phi->uses()) {
    //         if (!use.isReg()) continue;
    //         if (!TargetRegisterInfo::isVirtualRegister(use.getReg())) continue;
    //         outs() << "\tPHI: uses: "<< use << "\n";
    //     }
    //     outs() << "\n";
    //     outs() << "\n";
    // }

    ////////////////////////////////////////////////////////////////////////////////////
    // merge SSA LiveRanges for same base variables
    // use Phi node info already kept in this class
    // 
    outs() << "*** Total Number of live ranges: " << interferenceGraph.size() << "\n";
    for (MachineInstr* phi : getPhiNodes()) {
        // save live ranges for all virtual register operands 
        std::vector<LiveRange*> toMerge;
        for (unsigned int i = 0; i < phi->getNumOperands(); ++i) {
            if (!phi->getOperand(i).isReg()) continue;
            if (!TargetRegisterInfo::isVirtualRegister(phi->getOperand(i).getReg())) continue;

            LiveRange* LR = getLiveRangeFor(phi->getOperand(i).getReg());
            toMerge.push_back(LR);
        }

        // merge
        // remove merged LiveRange from interferenceGraph
        // update the interferenceGraphRetrieval map 
        // no need to update the MI liveness info structure because it is linked to 
        // reg numbers and not live ranges
        if (toMerge.size() > 1) {
            for (size_t i = 1; i < toMerge.size(); ++i) {
                outs() << "*** Merging ranges for reg #'s: ";
                for (unsigned regNum : toMerge[0]->virtRegs) {
                    outs() << TargetRegisterInfo::virtReg2Index(regNum) << " ";
                }
                outs() << " with reg #'s: ";
                for (unsigned regNum : toMerge[i]->virtRegs) {
                    outs() << TargetRegisterInfo::virtReg2Index(regNum) << " ";
                }
                outs() << "\n";

                toMerge[0]->MergeWith(toMerge[i]);
                interferenceGraph.erase(toMerge[i]);
                for (unsigned virtReg : toMerge[i]->virtRegs) {
                    interferenceGraphRetrieval[virtReg] = toMerge[0];
                }

                delete toMerge[i];
            }

        }
    }  // end live range (phi) merge loop
    outs() << "*** Total Number of live ranges: " << interferenceGraph.size() << "\n";
}

void X86RegisterAllocator::calcLivenessForDef(MachineInstr* def, MachineOperand& defOp, MachineRegisterInfo& regInfo) {
    // go up the CFG until we reach the def 
    // SSA means that def dominates all uses so we dont need to 
    // go over the whole cfg
    // for all such paths and all such instructions
    // mark the virtual register as live-in for those instructions

    std::set<MachineInstr*> liveInstrs;

    // share this set for all `uses` to prune search early
    // and avoid redundant work
    std::set<MachineBasicBlock*> VisitedBBlist;

    unsigned regNum = defOp.getReg();
    LiveRange* LR = getLiveRangeFor(regNum);

    // outs() << "\tWalking the CFG for all uses\n";
    for (MachineInstr& u : regInfo.use_instructions(regNum)) {
        MachineInstr* user = &u;
        MachineInstr* CurrMI = user;

        // Phi nodes are special. the use operands might not be dominated by their defs 
        // so propogate the liveness information 
        // across ONLY the corresponding parent branches
        if (user->isPHI()) {
            addPhiNode(user);

            // set CurrMI to be the last instruction of the phi's parent
            auto regOperandNum = user->findRegisterUseOperandIdx(regNum);
            auto MBBOperand = user->getOperand(regOperandNum + 1);
            auto MBB = MBBOperand.getMBB();
            CurrMI = &(MBB->instr_back());
        }

        std::queue<MachineBasicBlock*> BBlist;
        liveInstrs.insert(CurrMI);
        VisitedBBlist.insert(CurrMI->getParent());

        // liveness analysis starting for instrs in the current BB where we have the use
        while (MachineInstr* MI = CurrMI->getPrevNode()) {
            CurrMI = MI;
            if (MI == def) break;

            liveInstrs.insert(CurrMI);
        }

        // add all preds of current BB to worklist if not already visited
        if (CurrMI != def) {
            for (MachineBasicBlock* MBB : CurrMI->getParent()->predecessors()) {
                if (!VisitedBBlist.count(MBB)) BBlist.push(MBB);
            }
        }


        // liveness analysis for all other BBs
        while (!BBlist.empty()) {
            // get the last instr
            MachineBasicBlock* curBB = BBlist.front();
            BBlist.pop();

            if (VisitedBBlist.count(curBB)) {
                continue;
            }

            VisitedBBlist.insert(curBB);
            
            if (curBB->empty()) {
                continue;
            }

            CurrMI = &(curBB->instr_back());

            // reverse iterate
            for (MachineInstr* MI = CurrMI; MI != NULL; MI = CurrMI->getPrevNode()) {
                CurrMI = MI;
                if (MI == def) break;

                liveInstrs.insert(CurrMI);
            }

            // add all preds of current BB to worklist if not already visited
            if (CurrMI != def) {
                for (MachineBasicBlock* MBB : CurrMI->getParent()->predecessors()) {
                    if (!VisitedBBlist.count(MBB)) BBlist.push(MBB);
                }
            }
        }  // goto BB worklist's next element
    }

    // outs() << "\tUpdating Liveness data structures\n";

    // remove all phi's from the live set
    for (std::set<MachineInstr*>::iterator i = liveInstrs.begin(); i != liveInstrs.end();) {
        if ((*i)->isPHI()) {
            i = liveInstrs.erase(i);
        } else {
            ++i;
        }
    }

    // update the live range
    LR->addLiveInstrs(liveInstrs);

    // update the MI-to-virtRegs liveness map
    for (MachineInstr* MI : liveInstrs) {
        livenessInformation[MI].insert(regNum);
    }
}

} // end of anonymous namespace

char X86RegisterAllocator::ID = 0;

INITIALIZE_PASS(X86RegisterAllocator,
    X86_REGISTER_ALLOCATOR_PASS_NAME,
    X86_REGISTER_ALLOCATOR_PASS_DESC,
    true, // is CFG only?
    false  // is analysis?
)

FunctionPass *llvm::createX86RegisterAllocator() { 
    return new X86RegisterAllocator(); 
}
