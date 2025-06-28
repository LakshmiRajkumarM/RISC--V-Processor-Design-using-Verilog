`timescale 1ns / 1ps


// PROGRAM COUNTER
module Program_Counter(clk, rst, en, PC_next, PC);
    input clk, rst, en;
    input [31:0] PC_next;
    output reg [31:0] PC;

    always @(posedge clk or posedge rst) begin
        if (rst)
            PC <= 32'b0;
        else if (en)  
            PC <= PC_next;
    end
endmodule

// PC + 4
module PCPlus4(fromPC, ToPC);
    input [31:0] fromPC;
    output [31:0] ToPC;

    assign ToPC = fromPC + 4;
endmodule

// MUX PCNext
module muxPC_Next(PCPlus4, PCSrc, PCTarget, PCNext);
    input [31:0] PCPlus4, PCTarget;
    input PCSrc;
    output reg [31:0] PCNext;

    always @(*) begin
        case (PCSrc)
            1'b0: PCNext = PCPlus4;
            1'b1: PCNext = PCTarget;
            default: PCNext = PCPlus4;  
        endcase
    end
endmodule

// Instruction Memory
module Instruction_Memory(clk, rst, A, RD);
    input clk, rst;
    input [31:0] A;
    output [31:0] RD;  

    reg [31:0] Inst_mem[0:1023];  
    integer i;

    // Initialize memory on reset
    always @(posedge rst) begin
        if (rst) begin
            for (i = 0; i < 1024; i = i + 1)
                Inst_mem[i] <= 32'b0;
        end
    end

    assign RD = Inst_mem[A[11:2]];  
 endmodule

//************************************************************************************************************
// Fetch Register (IF/ID Pipeline Register)
module fetch_register(clk, en, clr, RD, PCF, PCPlus4F, InstrD, PCD, PCPlus4D);
    input clk, en, clr;
    input [31:0] RD, PCF, PCPlus4F;
    output reg [31:0] InstrD, PCD, PCPlus4D;

    always @(posedge clk) begin
        if (clr) begin
            InstrD <= 32'b0;
            PCD <= 32'b0;
            PCPlus4D <= 32'b0;
        end
        else if (en) begin  
            InstrD <= RD;
            PCD <= PCF;
            PCPlus4D <= PCPlus4F;
        end
    end
endmodule

//************************************************************************************************************
// Register File
module Register_file(clk, rst, we3, A1, A2, A3, RD1, RD2, WD3);
    input clk, rst, we3;
    input [4:0] A1, A2, A3;
    input [31:0] WD3;
    output [31:0] RD1, RD2;

    reg [31:0] Register[0:31];
    integer j;

    always @(negedge clk or posedge rst) begin
        if (rst) begin
            for (j = 0; j < 32; j = j + 1)
                Register[j] <= 32'b0;
        end else if (we3 && A3 != 0)
            Register[A3] <= WD3;
    end

    assign RD1 = Register[A1];
    assign RD2 = Register[A2];
endmodule


// Extend
module extend(input [31:0] instr, input [1:0] immsrc, output reg [31:0] immext);
    always @(*) begin
        case (immsrc)
            2'b00: immext = {{20{instr[31]}}, instr[31:20]};  // I-type
            2'b01: immext = {{20{instr[31]}}, instr[31:25], instr[11:7]};  // S-type
            2'b10: immext = {{20{instr[31]}}, instr[7], instr[30:25], instr[11:8], 1'b0};  // B-type
            2'b11: immext = {{12{instr[31]}}, instr[19:12], instr[20], instr[30:21], 1'b0};  // J-type
            default: immext = 32'b0;  
        endcase
    end
endmodule

//************************************************************************************************************
//Decode Register (ID/EX Pipeline Register)
module decode_register(clk, clr, RD1, RD2, PCD, Rs1D, Rs2D, RdD, ImmExtD, PCPlus4D, RegWriteD, ResultSrcD, MemWriteD, JumpD, BranchD,
                      ALUControlD, ALUSrcD, RD1E, RD2E, PCE, Rs1E, Rs2E, RdE, ImmExtE, PCPlus4E, RegWriteE, ResultSrcE, MemWriteE, JumpE, BranchE,
                      ALUControlE, ALUSrcE);
                        
    input clk, clr;
    input RegWriteD, MemWriteD, JumpD, BranchD, ALUSrcD;
    input [1:0] ResultSrcD;
    input [3:0] ALUControlD;
    input [31:0] RD1, RD2, PCD, PCPlus4D, ImmExtD;
    input [4:0] Rs1D, Rs2D, RdD;

    output reg [31:0] RD1E, RD2E, PCE, ImmExtE, PCPlus4E;                      
    output reg [4:0] Rs1E, Rs2E, RdE;
    output reg RegWriteE, MemWriteE, JumpE, BranchE, ALUSrcE;
    output reg [1:0] ResultSrcE;
    output reg [3:0] ALUControlE;

    always @(posedge clk) begin
        if (clr) begin
            RD1E <= 32'b0;
            RD2E <= 32'b0;
            PCE <= 32'b0;
            ImmExtE <= 32'b0;
            PCPlus4E <= 32'b0;
            Rs1E <= 5'b0;
            Rs2E <= 5'b0;
            RdE <= 5'b0;
            RegWriteE <= 1'b0;
            MemWriteE <= 1'b0;
            JumpE <= 1'b0;
            BranchE <= 1'b0;
            ALUSrcE <= 1'b0;
            ResultSrcE <= 2'b0;
            ALUControlE <= 4'b0;
        end
        else begin
            RD1E <= RD1;
            RD2E <= RD2;
            PCE <= PCD;
            ImmExtE <= ImmExtD;
            PCPlus4E <= PCPlus4D;
            Rs1E <= Rs1D;
            Rs2E <= Rs2D;
            RdE <= RdD;
            RegWriteE <= RegWriteD;
            MemWriteE <= MemWriteD;
            JumpE <= JumpD;
            BranchE <= BranchD;
            ALUSrcE <= ALUSrcD;
            ResultSrcE <= ResultSrcD;
            ALUControlE <= ALUControlD;
        end
    end
endmodule

//************************************************************************************************************
//FORWARD MUX 1
module forward_mux1(RD1E, ResultW, ALUResultM, ForwardAE, SrcAE);
    input [31:0] RD1E, ResultW, ALUResultM;
    input [1:0] ForwardAE;
    output reg [31:0] SrcAE;

    always @(*) begin
        case(ForwardAE)
            2'b00: SrcAE = RD1E;
            2'b01: SrcAE = ResultW;
            2'b10: SrcAE = ALUResultM;
            default: SrcAE = RD1E;
        endcase
    end
endmodule

//FORWARD MUX 2
module forward_mux2(RD2E, ResultW, ALUResultM, ForwardBE, SrcB_partial);
    input [31:0] RD2E, ResultW, ALUResultM;
    input [1:0] ForwardBE;
    output reg [31:0] SrcB_partial;

    always @(*) begin
        case(ForwardBE)
            2'b00: SrcB_partial = RD2E;
            2'b01: SrcB_partial = ResultW;
            2'b10: SrcB_partial = ALUResultM;
            default: SrcB_partial = RD2E;
        endcase
    end
endmodule

//HAZARD UNIT
module hazard_unit(
    input [4:0] Rs1E, Rs2E, RdM, RdW,
    input RegWriteM, RegWriteW, ResultSrcE0, PCSrcE,
    input [4:0] Rs1D, RdE, Rs2D,
    input rst,
    output reg [1:0] forwardAE, forwardBE,
    output reg stallF, stallD, flushD, flushE
);

wire lwstall;

// Forwarding logic for source A
always @(*) begin
    if (rst)
        forwardAE = 2'b00;
    else if (((Rs1E == RdM) && RegWriteM) && (Rs1E != 0))
        forwardAE = 2'b10;
    else if (((Rs1E == RdW) && RegWriteW) && (Rs1E != 0))
        forwardAE = 2'b01;
    else
        forwardAE = 2'b00;
end

// Forwarding logic for source B
always @(*) begin
    if (rst)
        forwardBE = 2'b00;
    else if (((Rs2E == RdM) && RegWriteM) && (Rs2E != 0))
        forwardBE = 2'b10;
    else if (((Rs2E == RdW) && RegWriteW) && (Rs2E != 0))
        forwardBE = 2'b01;
    else
        forwardBE = 2'b00;
end

// Load-use hazard detection
assign lwstall = ResultSrcE0 && ((Rs1D == RdE) || (Rs2D == RdE));

// Stall and flush control
always @(*) begin
    if (rst) begin
        stallF = 0;
        stallD = 0;
        flushD = 0;
        flushE = 0;
    end else begin
        stallF = lwstall;
        stallD = lwstall;
        flushD = PCSrcE;
        flushE = PCSrcE || lwstall;
    end
end

endmodule

// MUX
module mux2_1(SrcB_partial, ImmExt, ALUSrc, SrcB);  // Corrected: parameter order
    input [31:0] SrcB_partial, ImmExt;
    input ALUSrc;
    output [31:0] SrcB;

    assign SrcB = ALUSrc ? ImmExt : SrcB_partial;
endmodule

// ALU
module ALU (
    input  [31:0] srcA,
    input  [31:0] srcB,
    input  [3:0]  ALUControl,       
    output reg [31:0] ALUResult,
    output Zero
);
    always @(*) begin
        case (ALUControl)
            4'b0000: ALUResult = srcA + srcB;                                   // ADD
            4'b0001: ALUResult = srcA - srcB;                                   // SUB
            4'b0010: ALUResult = srcA & srcB;                                   // AND
            4'b0011: ALUResult = srcA | srcB;                                   // OR
            4'b0101: ALUResult = ($signed(srcA) < $signed(srcB)) ? 32'b1 : 32'b0; // SLT

            // M-extension Instructions
            4'b1000: ALUResult = $signed(srcA) * $signed(srcB);                 // MUL
            4'b1100: ALUResult = (srcB != 0) ? $signed(srcA) / $signed(srcB) : 32'hFFFFFFFF; // DIV
            4'b1110: ALUResult = (srcB != 0) ? $signed(srcA) % $signed(srcB) : srcA;         // REM

            default: ALUResult = 32'b0;
        endcase
    end

    assign Zero = (ALUResult == 32'b0);
endmodule

// PC Target
module PCTarget(PC, ImmExt, PCTargetOut);
    input [31:0] PC, ImmExt;
    output [31:0] PCTargetOut;

    assign PCTargetOut = PC + ImmExt;
endmodule

//************************************************************************************************************
//Execute Register (EX/MEM Pipeline Register)
module execute_register(clk, ALUResultE, WriteDataE, RdE, PCPlus4E, RegWriteE, ResultSrcE, MemWriteE, ALUResultM,
                       WriteDataM, RdM, PCPlus4M, RegWriteM, ResultSrcM, MemWriteM);

    input clk;
    input RegWriteE, MemWriteE;
    input [1:0] ResultSrcE;
    input [31:0] ALUResultE, WriteDataE, PCPlus4E;
    input [4:0] RdE;
    output reg [4:0] RdM;
    output reg [31:0] ALUResultM, WriteDataM, PCPlus4M;
    output reg RegWriteM, MemWriteM;
    output reg [1:0] ResultSrcM;

    always @(posedge clk) begin
        ALUResultM <= ALUResultE;
        WriteDataM <= WriteDataE;
        PCPlus4M <= PCPlus4E;
        RdM <= RdE;
        RegWriteM <= RegWriteE;
        MemWriteM <= MemWriteE;
        ResultSrcM <= ResultSrcE;
    end
endmodule

//************************************************************************************************************
// Data Memory
module Data_Memory(clk, rst, WE, A, WD, RD);
    input clk, rst, WE;
    input [31:0] A, WD;
    output [31:0] RD;  

    reg [31:0] Reg[0:1023];  
    integer k;

    always @(posedge clk or posedge rst) begin
        if (rst) begin
            for (k = 0; k < 1024; k = k + 1)
                Reg[k] <= 32'b0;
        end else if (WE)
            Reg[A[11:2]] <= WD;  
    end

    assign RD = Reg[A[11:2]];
endmodule

// Main Decoder
module maindec(
    input [6:0] op,
    output reg [1:0] ResultSrc,
    output reg MemWrite,
    output reg Branch,
    output reg ALUSrc,
    output reg RegWrite,
    output reg Jump,
    output reg [1:0] ImmSrc,
    output reg [1:0] ALUOp
);

    always @(*) begin
        case (op)
            7'b0110011: begin // R-type
                RegWrite = 1; ImmSrc = 2'bxx; ALUSrc = 0; MemWrite = 0;
                ResultSrc = 2'b00; Branch = 0; ALUOp = 2'b10; Jump = 0;
            end
            7'b0010011: begin // I-type (arithmetic)
                RegWrite = 1; ImmSrc = 2'b00; ALUSrc = 1; MemWrite = 0;
                ResultSrc = 2'b00; Branch = 0; ALUOp = 2'b10; Jump = 0;
            end
            7'b0000011: begin // Load
                RegWrite = 1; ImmSrc = 2'b00; ALUSrc = 1; MemWrite = 0;
                ResultSrc = 2'b01; Branch = 0; ALUOp = 2'b00; Jump = 0;
            end
            7'b0100011: begin // Store
                RegWrite = 0; ImmSrc = 2'b01; ALUSrc = 1; MemWrite = 1;
                ResultSrc = 2'bxx; Branch = 0; ALUOp = 2'b00; Jump = 0;
            end
            7'b1100011: begin // Branch
                RegWrite = 0; ImmSrc = 2'b10; ALUSrc = 0; MemWrite = 0;
                ResultSrc = 2'bxx; Branch = 1; ALUOp = 2'b01; Jump = 0;
            end
            7'b1101111: begin // JAL
                RegWrite = 1; ImmSrc = 2'b11; ALUSrc = 0; MemWrite = 0;
                ResultSrc = 2'b10; Branch = 0; ALUOp = 2'b00; Jump = 1;
            end
            default: begin
                RegWrite = 0; ImmSrc = 2'b00; ALUSrc = 0; MemWrite = 0;
                ResultSrc = 2'b00; Branch = 0; ALUOp = 2'b00; Jump = 0;
            end
        endcase
    end
endmodule

// ALU Decoder
module aludec (
    input  wire [1:0] ALUOp,
    input  wire [2:0] funct3,
    input  wire       opb5,       
    input  wire       funct7b5,funct7b1,   
    output reg  [3:0] ALUControl 
);
    always @(*) begin
        case (ALUOp)
            2'b00: ALUControl = 4'b0000; 
            2'b01: ALUControl = 4'b0001; 
            2'b10: begin // R-type or I-type arithmetic
                if (opb5 && funct7b1) begin
                    // R-type M-extension instructions
                    case (funct3)
                        3'b000: ALUControl = 4'b1000; // MUL
                        3'b100: ALUControl = 4'b1100; // DIV
                        3'b110: ALUControl = 4'b1110; // REM
                        default: ALUControl = 4'b1000;
                    endcase
                end else begin
                    // Regular arithmetic
                    case (funct3)
                        3'b000: ALUControl = (opb5 && funct7b5) ? 4'b0001 : 4'b0000; // SUB or ADD
                        3'b010: ALUControl = 4'b0101; // SLT
                        3'b110: ALUControl = 4'b0011; // OR
                        3'b111: ALUControl = 4'b0010; // AND
                        default: ALUControl = 4'b0000;
                    endcase
                end
            end
            default: ALUControl = 4'b0000;
        endcase
    end
endmodule

// Top-Level Controller
module controller(
    input [6:0] op,
    input [2:0] funct3,
  input funct7,funct7_1,
    input Zero,
    output [1:0] ResultSrc,
    output MemWrite,
    output ALUSrc,
    output RegWrite,
    output Jump,
    output [1:0] ImmSrc,
    output [3:0] ALUControl,
    output Branch
);
    wire [1:0] ALUOp;
  	wire funct7b5 = funct7;
    wire funct7b1 = funct7_1;


    maindec md(
        .op(op), .ResultSrc(ResultSrc), .MemWrite(MemWrite), .Branch(Branch),
        .ALUSrc(ALUSrc), .RegWrite(RegWrite), .Jump(Jump), .ImmSrc(ImmSrc), .ALUOp(ALUOp)
    );

    aludec ad(
      .opb5(op[5]), .funct3(funct3),.funct7b5(funct7b5),.funct7b1(funct7b1), .ALUOp(ALUOp), .ALUControl(ALUControl)
    );
endmodule

//************************************************************************************************************
//Memory Register (MEM/WB Pipeline Register)
module memory_register(clk, RD, RdM, PCPlus4M, RegWriteM, ResultSrcM, ALUResultM, ReadData, RdW, PCPlus4W, RegWriteW, ResultSrcW, ALUResultW);
    input clk;
    input [31:0] RD, PCPlus4M, ALUResultM;
    input [4:0] RdM;
    input RegWriteM;
    input [1:0] ResultSrcM;
    output reg [31:0] ReadData, PCPlus4W, ALUResultW;
    output reg [4:0] RdW;
    output reg RegWriteW;
    output reg [1:0] ResultSrcW;

    always @(posedge clk) begin
        ReadData <= RD;
        PCPlus4W <= PCPlus4M;  
        ALUResultW <= ALUResultM;
        RdW <= RdM;
        RegWriteW <= RegWriteM;
        ResultSrcW <= ResultSrcM;
    end
endmodule

//************************************************************************************************************
// Final MUX
module mux_last(ALUResult, ReadData, PC_plus4, ResultSrc, Result);
    input [31:0] ALUResult, ReadData, PC_plus4;
    input [1:0] ResultSrc;
    output reg [31:0] Result;

    always @(*) begin
        case (ResultSrc)
            2'b00: Result = ALUResult;
            2'b01: Result = ReadData;
            2'b10: Result = PC_plus4;
            default: Result = ALUResult;  // Corrected: default to ALUResult
        endcase
    end
endmodule

// Top Module
module PIPELINE_PROCESSOR(
    input clk, rst,
    output [31:0] PC,
    output [31:0] Instr,
    output MemWriteM, stallF, 
    output [31:0] ALUResultM, WriteDataM, ResultW, RD1, SrcAE, SrcBE,
    output [3:0] ALUControlE
);
    
    wire [31:0] PC_next, PCPlus4F, PCTargetE, ImmExtD;
    wire [31:0] RD2, RD, PCD, PCPlus4D, RDD;
    wire [31:0] RD1E, RD2E, PCE, ImmExtE, PCPlus4E;
    wire [31:0] SrcB_partial, ALUResultE, WriteDataE;
    wire [31:0] PCPlus4M, ReadDataW, PCPlus4W, ALUResultW;
    wire [4:0] Rs1D, Rs2D, RdD, Rs1E, Rs2E, RdE, RdM, RdW;
    wire [1:0] ResultSrcD, ImmSrcD, ForwardAE, ForwardBE;
    wire [1:0] ResultSrcE, ResultSrcM, ResultSrcW;
    wire [3:0] ALUControlD;
    wire RegWriteD, MemWriteD, JumpD, BranchD, ALUSrcD;
    wire RegWriteE, MemWriteE, JumpE, BranchE, ALUSrcE;
    wire RegWriteM, RegWriteW;
    wire PCSrcE, ZeroE, ResultSrcE0;
    wire stallD, flushD, flushE;

    assign Rs1D = Instr[19:15];
    assign Rs2D = Instr[24:20];
    assign RdD = Instr[11:7];
    assign ResultSrcE0 = ResultSrcE[0];
    assign WriteDataE = SrcB_partial; 
    assign PCSrcE = ((BranchE & ZeroE) | JumpE);

    Program_Counter m1(clk, rst, ~stallF, PC_next, PC);  
    PCPlus4 m2(PC, PCPlus4F);
    muxPC_Next m3(PCPlus4F, PCSrcE, PCTargetE, PC_next);
    Instruction_Memory im1(clk, rst, PC, RD);
  
    fetch_register m4(clk, ~stallD, flushD, RD, PC, PCPlus4F, Instr, PCD, PCPlus4D);  // Corrected: enable when not stalled
    Register_file m5(clk, rst, RegWriteW, Rs1D, Rs2D, RdW, RD1, RD2, ResultW);
  controller m6(Instr[6:0], Instr[14:12], Instr[30],Instr[25], ZeroE, ResultSrcD, MemWriteD, ALUSrcD, RegWriteD, JumpD, ImmSrcD, ALUControlD, BranchD);
    extend m7(Instr, ImmSrcD, ImmExtD);
    decode_register m8(clk, flushE, RD1, RD2, PCD, Rs1D, Rs2D, RdD, ImmExtD, PCPlus4D, RegWriteD, ResultSrcD, MemWriteD, JumpD, BranchD,
                      ALUControlD, ALUSrcD, RD1E, RD2E, PCE, Rs1E, Rs2E, RdE, ImmExtE, PCPlus4E, RegWriteE, ResultSrcE, MemWriteE, JumpE, BranchE,
                      ALUControlE, ALUSrcE);
    forward_mux1 m9(RD1E, ResultW, ALUResultM, ForwardAE, SrcAE);
    forward_mux2 m10(RD2E, ResultW, ALUResultM, ForwardBE, SrcB_partial);
    mux2_1 m11(SrcB_partial, ImmExtE, ALUSrcE, SrcBE);
    PCTarget m12(PCE, ImmExtE, PCTargetE);
    ALU m13(SrcAE, SrcBE, ALUControlE, ALUResultE, ZeroE);
    execute_register m14(clk, ALUResultE, WriteDataE, RdE, PCPlus4E, RegWriteE, ResultSrcE, MemWriteE, ALUResultM,
                        WriteDataM, RdM, PCPlus4M, RegWriteM, ResultSrcM, MemWriteM);
    Data_Memory m15(clk, rst, MemWriteM, ALUResultM, WriteDataM, RDD);
    memory_register m16(clk, RDD, RdM, PCPlus4M, RegWriteM, ResultSrcM, ALUResultM, ReadDataW, RdW, PCPlus4W, RegWriteW, ResultSrcW, ALUResultW);
    mux_last m17(ALUResultW, ReadDataW, PCPlus4W, ResultSrcW, ResultW);
    hazard_unit m18(Rs1E, Rs2E, RdM, RdW, RegWriteM, RegWriteW, ResultSrcE0, PCSrcE,
                   Rs1D, RdE, Rs2D, rst, ForwardAE, ForwardBE, stallF, stallD, flushD, flushE);

endmodule
