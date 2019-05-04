#include "X86.h"
#include "X86InstrInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstrBuilder.h"
#include "llvm/CodeGen/MachineRegisterInfo.h"

using namespace llvm;

#define X86_MACHINEINSTR_PRINTER_PASS_DESC "Dummy X86 machineinstr printer pass"
#define X86_MACHINEINSTR_PRINTER_PASS_NAME "x86-machineinstr-printer"

namespace {

class X86MachineInstrPrinter : public MachineFunctionPass {
public:
    static char ID;

    X86MachineInstrPrinter() : MachineFunctionPass(ID) {
        initializeX86MachineInstrPrinterPass(*PassRegistry::getPassRegistry());
    }

    bool runOnMachineFunction(MachineFunction &MF) override;

    StringRef getPassName() const override { return X86_MACHINEINSTR_PRINTER_PASS_NAME; }
};

bool X86MachineInstrPrinter::runOnMachineFunction(MachineFunction &MF) {

    for (auto &MBB : MF) {

        outs() << "Contents of MachineBasicBlock:\n";
        outs() << MBB << "\n";

        // const BasicBlock *BB = MBB.getBasicBlock();
        // outs() << "\n";
        // outs() << "Contents of BasicBlock corresponding to MachineBasicBlock:\n";
        // outs() << *BB << "\n";
        // outs() << "\n\n";
    }
        outs() << "\n\n";

    // MachineRegisterInfo& regInfo = MF.getRegInfo();//(&MF);
    // outs() << "Printing Instructions of MachineFunction:\n";
    // for (auto &MBB : MF) {
    //     for (auto& MI : MBB) {
    //         outs() << MI;
    //         for (MachineOperand& def : MI.defs()) {
    //             if (!def.isReg()) continue;
    //             if (!TargetRegisterInfo::isVirtualRegister(def.getReg())) continue;

    //             unsigned regNum = def.getReg();
    //             outs() << "\tdef: "<< def 
    //             << " Register Number: " << (unsigned)TargetRegisterInfo::virtReg2Index(regNum) 
    //             << "\n";
    //             for (auto& users : regInfo.use_instructions(def.getReg())) {
    //                 outs() << "\t\tuse: "<< users << "\n";
    //             }
    //             for (auto& users : regInfo.def_instructions(def.getReg())) {
    //                 outs() << "\t\tdef: "<< users << "\n";
    //             }                
    //         }
    //         // for (unsigned i = 0, e = MRI->getNumVirtRegs(); i != e; ++i) {
    //         //   // reg ID
    //         //   unsigned Reg = TargetRegisterInfo::index2VirtReg(i);

    //         // for (auto& uses : MI.uses()) {
    //         //     outs() << "\tuses: "<< uses << "\n";
    //         //     // for (auto& users : uses.users()) {
    //         //     //     outs() << "\tusers: "<< users << "\n";
    //         //     // }
    //         // }
    //         outs() << "\n";
    //     }
    //     outs() << "\n\n";
    // }    
    return false;
}

} // end of anonymous namespace

char X86MachineInstrPrinter::ID = 0;

INITIALIZE_PASS(X86MachineInstrPrinter,
    X86_MACHINEINSTR_PRINTER_PASS_NAME,
    X86_MACHINEINSTR_PRINTER_PASS_DESC,
    true, // is CFG only?
    true  // is analysis?
)

FunctionPass *llvm::createX86MachineInstrPrinter() { 
    return new X86MachineInstrPrinter(); 
}
