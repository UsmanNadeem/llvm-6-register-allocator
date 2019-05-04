#include "X86.h"
#include "X86InstrInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include "llvm/CodeGen/MachineRegisterInfo.h"

#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <queue>

using namespace llvm;

#define X86_REGISTER_ALLOCATOR_PASS_DESC "A basic register allocator for X86"
#define X86_REGISTER_ALLOCATOR_PASS_NAME "x86-reg-alloc"

namespace {

class LiveRange {
public:
    // std::vector<MachineBasicBlock*> BBs;
    // list of instructions where this equivalence class is live-in
    std::set<MachineInstr*> liveInstrs;

    // regs in one  equivalence class
    std::set<unsigned> virtRegs;
    
    // the physical register allocated to this LiveRange
    unsigned physicalReg;
    bool isAllocated;

    LiveRange() {isAllocated = false;}
    LiveRange(unsigned virtReg) {
        isAllocated = false;
        virtRegs.insert(virtReg);
    }
    ~LiveRange() {}


    bool contains (unsigned virtReg) {
        return virtRegs.count(virtReg);
    }
    
    bool contains (MachineInstr* MI) {
        return liveInstrs.count(MI);
    }    
    void MergeWith (LiveRange* LR) {
        liveInstrs.insert(LR->liveInstrs.begin(),LR->liveInstrs.end());
        virtRegs.insert(LR->virtRegs.begin(),LR->virtRegs.end());
    }

    bool interferesWith (LiveRange* LR) {
        std::vector<MachineInstr*> unionSet;

        std::set_union (liveInstrs.begin(), liveInstrs.end(), 
            LR->liveInstrs.begin(), LR->liveInstrs.end(), unionSet.begin());

        return (unionSet.size() > 0);
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

    X86RegisterAllocator() : MachineFunctionPass(ID) {
        initializeX86RegisterAllocatorPass(*PassRegistry::getPassRegistry());
    }

    bool runOnMachineFunction(MachineFunction &MF) override;

    StringRef getPassName() const override { 
        return X86_REGISTER_ALLOCATOR_PASS_NAME; 
    }

    void calcGlobalLivenessInfo(MachineFunction &MF, MachineRegisterInfo& regInfo);

    // .. todo merge ranges ..
    // might want to change to a vector
    std::set<LiveRange*> interferenceGraph;
    std::map<unsigned, LiveRange*> interferenceGraphRetrieval;
    std::map<MachineInstr*, std::set<unsigned>> livenessInformation;

    // list of phi nodes that we need to remove after reg allocation
    std::set<MachineInstr*> phiNodes;

    std::set<MachineInstr*> getPhiNodes() {
        return phiNodes;
    }
    
    void addPhiNode (MachineInstr* MI) {
        phiNodes.insert(MI);
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

    void calcLivenessForDef(MachineInstr* def, MachineOperand& defOp, MachineRegisterInfo& regInfo) {
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

        outs() << "\tWalking the CFG for all uses\n";
        for (MachineInstr& u : regInfo.use_instructions(regNum)) {
            MachineInstr* user = &u;
            MachineInstr* CurrMI = user;

            // Phi nodes are special. the use operands might not be dominated by their defs 
            // so propogate the liveness information 
            // across ONLY the corresponding parent branches
            if (user->isPHI()) {
                addPhiNode(user);

                // set CurrMI to be the last instruction of the phi's parent
                CurrMI = &(user->getOperand(user->findRegisterUseOperandIdx(regNum) + 1).getMBB()->instr_back());
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

        outs() << "\tUpdating Liveness data structures\n";

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

};

bool X86RegisterAllocator::runOnMachineFunction(MachineFunction &MF) {
    // cleanup from previous function run...if any
    for (LiveRange* LR : interferenceGraph) {
        delete LR;
    }
    interferenceGraph.clear();
    interferenceGraphRetrieval.clear();
    livenessInformation.clear();
    phiNodes.clear();

    MachineRegisterInfo& regInfo = MF.getRegInfo();
    calcGlobalLivenessInfo(MF, regInfo);

    return false;
}

void X86RegisterAllocator::calcGlobalLivenessInfo(MachineFunction &MF, MachineRegisterInfo& regInfo) {
// .. todo look at already live physicalRegs ..
    
    // liveness analysis for every def/virtual SSA Reg 
    for (auto &MBB : MF) {
        for (auto& MI : MBB) {
            for (MachineOperand& def : MI.defs()) {
                if (!def.isReg()) continue;
                if (!TargetRegisterInfo::isVirtualRegister(def.getReg())) continue;
                    
                // perform Liveness analysis and 
                // update the LiveRange for the virtReg, and
                // update the MI to virtRegs liveness map
                unsigned regNum = def.getReg();
                outs() << "*** calcLivenessForDef reg #" << TargetRegisterInfo::virtReg2Index(regNum) << "\n";

                calcLivenessForDef(&MI, def, regInfo);
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

} // end of anonymous namespace

char X86RegisterAllocator::ID = 0;

INITIALIZE_PASS(X86RegisterAllocator,
    X86_REGISTER_ALLOCATOR_PASS_NAME,
    X86_REGISTER_ALLOCATOR_PASS_DESC,
    true, // is CFG only?
    true  // is analysis?
)

FunctionPass *llvm::createX86RegisterAllocator() { 
    return new X86RegisterAllocator(); 
}
