`timescale 1ns / 1ps


// PROGRAM COUNTER
module Program_Counter(clk, rst, PC_next, PC);
    input clk, rst;
    input [31:0] PC_next;
    output reg [31:0] PC;

    always @(posedge clk or posedge rst) begin
        if (rst)
            PC <= 32'b0;
        else
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
            default: PCNext = 32'b0;
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

    
    assign RD = Inst_mem[A[11:2]];  // Word-aligned addressing
endmodule

// Register File
module Register_file(clk, rst, we3, A1, A2, A3, RD1, RD2, WD3);
    input clk, rst, we3;
    input [4:0] A1, A2, A3;
    input [31:0] WD3;
    output [31:0] RD1, RD2;

    reg [31:0] Register[0:31];
    integer j;

    always @(posedge clk or posedge rst) begin
        if (rst) begin
            for (j = 0; j < 32; j = j + 1)
                Register[j] <= 32'b0;
        end else if (we3 && A3 != 0)
            Register[A3] <= WD3;
    end

    assign RD1 = Register[A1];
    assign RD2 = Register[A2];
endmodule


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

    // Combinational read
    assign RD = Reg[A[11:2]];
endmodule

// Extend
module extend(input [31:0] instr, input [1:0] immsrc, output reg [31:0] immext);
    always @(*) begin
        case (immsrc)
            2'b00: immext = {{20{instr[31]}}, instr[31:20]};
            2'b01: immext = {{20{instr[31]}}, instr[31:25], instr[11:7]};
            2'b10: immext = {{20{instr[31]}}, instr[7], instr[30:25], instr[11:8], 1'b0};
            2'b11: immext = {{12{instr[31]}}, instr[19:12], instr[20], instr[30:21], 1'b0};
            default: immext = 32'bx;
        endcase
    end
endmodule

// MUX
module mux2_1(RD2, ImmExt, ALUSrc, SrcB);
    input [31:0] RD2, ImmExt;
    input ALUSrc;
    output [31:0] SrcB;

    assign SrcB = ALUSrc ? ImmExt : RD2;
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
            4'b0000: ALUResult = srcA + srcB;                                   
            4'b0001: ALUResult = srcA - srcB;                                    
            4'b0010: ALUResult = srcA & srcB;                                    
            4'b0011: ALUResult = srcA | srcB;                                   
            4'b0101: ALUResult = ($signed(srcA) < $signed(srcB)) ? 32'b1 : 32'b0;

            // M-extension Instructions
            4'b1000: ALUResult = (($signed(srcA)) * ($signed(srcB)));                 // MUL
            4'b1100: ALUResult = (srcB != 0) ? $signed(srcA) / $signed(srcB) : 32'hFFFFFFFF; // DIV
            4'b1110: ALUResult = (srcB != 0) ? $signed(srcA) % $signed(srcB) : srcA;         // REM

            default: ALUResult = 32'b0;
        endcase
    end

    assign Zero = (ALUResult == 32'b0);
endmodule

// PC Target
module PCTarget(PC, ImmExt, PCTarget);
    input [31:0] PC, ImmExt;
    output [31:0] PCTarget;

    assign PCTarget = PC + ImmExt;
endmodule

// Main Decoder
module maindec(
    input [6:0] op,
    output [1:0] ResultSrc,
    output MemWrite,
    output Branch,
    output ALUSrc,
    output RegWrite,
    output Jump,
    output [1:0] ImmSrc,
    output [1:0] ALUOp
);
    reg [10:0] controls;
    assign {RegWrite, ImmSrc, ALUSrc, MemWrite, ResultSrc, Branch, ALUOp, Jump} = controls;

    always @(*) begin
        case (op)
            7'b0110011: controls = 11'b1_xx_0_0_00_0_10_0;
            7'b0010011: controls = 11'b1_00_1_0_00_0_10_0;
            7'b0000011: controls = 11'b1_00_1_0_01_0_00_0;
            7'b0100011: controls = 11'b0_01_1_1_00_0_00_0;
            7'b1100011: controls = 11'b0_10_0_0_00_1_01_0;
            7'b1101111: controls = 11'b1_11_0_0_10_0_00_1;
            default:    controls = 11'b0_00_0_0_00_0_00_0;
        endcase
    end
endmodule

// ALU Decoder
module aludec (
    input  wire [1:0] ALUOp,
    input  wire [2:0] funct3,
    input  wire       opb5,       
    input  wire       funct7b5,   
    output reg  [3:0] ALUControl 
);
    always @(*) begin
        case (ALUOp)
            2'b00: ALUControl = 4'b0000; 
            2'b01: ALUControl = 4'b0001; 
            2'b10: begin // R-type or I-type arithmetic
                if (opb5 && funct7b5) begin
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
    input funct7b5,
    input Zero,
    output [1:0] ResultSrc,
    output MemWrite,
    output PCSrc,
    output ALUSrc,
    output RegWrite,
    output Jump,
    output [1:0] ImmSrc,
    output [3:0] ALUControl
);
    wire [1:0] ALUOp;
    wire Branch;

    maindec md(
        .op(op), .ResultSrc(ResultSrc), .MemWrite(MemWrite), .Branch(Branch),
        .ALUSrc(ALUSrc), .RegWrite(RegWrite), .Jump(Jump), .ImmSrc(ImmSrc), .ALUOp(ALUOp)
    );

    aludec ad(
        .opb5(op[5]), .funct3(funct3), .funct7b5(funct7b5), .ALUOp(ALUOp), .ALUControl(ALUControl)
    );

    assign PCSrc = (Branch & Zero) | Jump;
endmodule

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
            default: Result = 32'b0;
        endcase
    end
endmodule

// Top Module
module rv32i(input clk, rst,
    output [31:0] PC,
    input [31:0] Instr,
    output MemWrite,
    output [31:0] ALUResult, WriteData, Result,RD1,SrcB,output [3:0]ALUControl
);
    wire [31:0] PC_next, PC_plus4, PCTarget, immext;
    wire [31:0]RD2, RD;
    wire [1:0] ResultSrc, ImmSrc;
    wire ALUSrc, RegWrite, Jump, PCSrc, zero;
    //wire [3:0] ALUControl;
    
    //output [31:0] RD1,SrcB;

    Program_Counter m1(clk, rst, PC_next, PC);
    PCPlus4 m2(PC, PC_plus4);
    muxPC_Next m3(PC_plus4, PCSrc, PCTarget, PC_next);
    Register_file m5(clk, rst, RegWrite, Instr[19:15], Instr[24:20], Instr[11:7], RD1, RD2, Result);
    extend m6(Instr, ImmSrc, immext);
    PCTarget m7(PC, immext, PCTarget);
    mux2_1 m8(RD2, immext, ALUSrc, SrcB);
    ALU m9(RD1, SrcB, ALUControl, ALUResult, zero);
    Data_Memory m10(clk, rst, MemWrite, ALUResult, RD2, RD);
    controller m11(Instr[6:0], Instr[14:12], Instr[25], zero, ResultSrc, MemWrite, PCSrc, ALUSrc, RegWrite, Jump, ImmSrc, ALUControl);
    mux_last m12(ALUResult, RD, PC_plus4, ResultSrc, Result);

    assign WriteData = RD2;
endmodule
