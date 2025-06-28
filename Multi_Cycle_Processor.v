// Code your design here
`timescale 1ns / 1ps


// PROGRAM COUNTER
module Program_Counter(clk, rst, PCWrite, PC_next, PC);
    input clk, rst, PCWrite;
    input [31:0] PC_next;
    output reg [31:0] PC;

    always @(posedge clk or posedge rst) begin
        if (rst)
            PC <= 32'b0;
        else if(PCWrite)
            PC <= PC_next;
    end
endmodule

// Instruction Register
module instruction_register(clk, rst, PC, IRWrite, RD, Instr, OldPC);
    input clk, rst, IRWrite;
    input [31:0] RD, PC;
    output reg [31:0] Instr, OldPC;

    always @(posedge clk or posedge rst) begin
        if(rst) begin
            Instr <= 32'b0;
            OldPC <= 32'b0;
        end
        else if(IRWrite) begin
            Instr <= RD;
            OldPC <= PC;
        end
    end
endmodule

// Address MUX
module mux2_1(PC, Result, AdrSrc, Adr);
    input [31:0] PC, Result;
    input AdrSrc;
    output reg [31:0] Adr;

    always @(*) begin
        case(AdrSrc)
            1'b0: Adr = PC;
            1'b1: Adr = Result;
        endcase
    end
endmodule

// Instruction/Data Memory (Combined)
module Instruction_Memory(clk, rst, MemWrite, A, WD, RD);
    input clk, rst, MemWrite;
    input [31:0] A;
    input [31:0] WD;
    output reg [31:0] RD;

    reg [31:0] Memory[0:1023]; // Increased memory size
    integer i;

    // Initialize memory with some test instructions
    initial begin
        
        for (i = 4; i < 1024; i = i + 1)
            Memory[i] = 32'h00000013; // Fill rest with NOPs
    end

    always @(posedge clk or posedge rst) begin
        if (rst) begin
            // Memory initialization handled in initial block
        end 
        else if(MemWrite) begin
            Memory[A[11:2]] <= WD; // Word-aligned access
        end
    end
    
    // Combinational read
    always @(*) begin
        RD = Memory[A[11:2]]; // Word-aligned access
    end
endmodule

// Data Register
module Data_Register(clk, rst, ReadData, Data);
    input clk, rst;
    input [31:0] ReadData;
    output reg [31:0] Data;

    always @(posedge clk or posedge rst) begin
        if(rst)
            Data <= 32'b0;
        else
            Data <= ReadData;
    end
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

    assign RD1 = (A1 == 0) ? 32'b0 : Register[A1];
    assign RD2 = (A2 == 0) ? 32'b0 : Register[A2];
endmodule

// A and WriteData Register
module rd1_register(clk, rst, RD1, RD2, A, WriteData);
    input clk, rst;
    input [31:0] RD1, RD2;
    output reg [31:0] A, WriteData;

    always @(posedge clk or posedge rst) begin
        if(rst) begin
            A <= 32'b0;
            WriteData <= 32'b0;
        end
        else begin
            A <= RD1;
            WriteData <= RD2;
        end
    end
endmodule

// ALU Source A MUX
module mux1(PC, OldPC, A, ALUSrcA, SrcA);
    input [31:0] PC, OldPC, A;
    input [1:0] ALUSrcA;
    output reg [31:0] SrcA;

    always @(*) begin
        case (ALUSrcA)
            2'b00: SrcA = PC;
            2'b01: SrcA = OldPC;
            2'b10: SrcA = A;
            default: SrcA = A;
        endcase
    end
endmodule

// ALU Source B MUX
module mux2(WriteData, ImmExt, ALUSrcB, SrcB);
    input [31:0] WriteData, ImmExt;
    input [1:0] ALUSrcB;
    output reg [31:0] SrcB;

    always @(*) begin
        case (ALUSrcB)
            2'b00: SrcB = WriteData;
            2'b01: SrcB = ImmExt;
            2'b10: SrcB = 32'd4;
            default: SrcB = 32'b0; 
        endcase
    end
endmodule
 
// Result MUX
module mux3(ALUOut, Data, ALUResult, ResultSrc, Result);
    input [31:0] ALUOut, Data, ALUResult;
    input [1:0] ResultSrc;
    output reg [31:0] Result;

    always @(*) begin
        case(ResultSrc)
            2'b00: Result = ALUOut;
            2'b01: Result = Data;
            2'b10: Result = ALUResult;
            default: Result = ALUOut;
        endcase
    end
endmodule

// Sign Extend
module extend(input [31:0] instr, input [1:0] immsrc, output reg [31:0] immext);
    always @(*) begin
        case (immsrc)
            2'b00: immext = {{20{instr[31]}}, instr[31:20]};                    // I-type
            2'b01: immext = {{20{instr[31]}}, instr[31:25], instr[11:7]};       // S-type
            2'b10: immext = {{20{instr[31]}}, instr[7], instr[30:25], instr[11:8], 1'b0}; // B-type
            2'b11: immext = {{12{instr[31]}}, instr[19:12], instr[20], instr[30:21], 1'b0}; // J-type
            default: immext = 32'b0;
        endcase
    end
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

// ALU Output Register
module ALU_register(clk, rst, ALUResult, ALUOut);
    input clk, rst;
    input [31:0] ALUResult;
    output reg [31:0] ALUOut;

    always @(posedge clk or posedge rst) begin
        if(rst)
            ALUOut <= 32'b0;
        else
            ALUOut <= ALUResult;
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


// Instruction Decoder
module Instr_Decoder(
    input  wire [6:0] op,
    output reg  [1:0] ImmSrc
);
    always @(*) begin
        case (op)
            7'b0000011: ImmSrc = 2'b00; // lw - I-type
            7'b0100011: ImmSrc = 2'b01; // sw - S-type
            7'b0110011: ImmSrc = 2'b00; // R-type 
            7'b1100011: ImmSrc = 2'b10; // beq - B-type
            7'b0010011: ImmSrc = 2'b00; // addi - I-type
            7'b1101111: ImmSrc = 2'b11; // jal - J-type
            default:    ImmSrc = 2'b00;
        endcase
    end
endmodule

// Main FSM
module main_fsm (
    input clk, rst, zero,
    input [6:0] op,
    output reg Branch, PCUpdate, PCWrite, RegWrite, MemWrite, IRWrite, AdrSrc,
    output reg [1:0] ResultSrc, ALUSrcA, ALUSrcB, ALUOp
);

    parameter S0_FETCH    = 4'd0,
              S1_DECODE   = 4'd1,
              S2_MEMADR   = 4'd2,
              S3_MEMREAD  = 4'd3,
              S4_MEMWB    = 4'd4,
              S5_MEMWRITE = 4'd5,
              S6_EXECUTER = 4'd6,
              S7_ALUWB    = 4'd7,
              S8_EXECUTEI = 4'd8,
              S9_JAL      = 4'd9,
              S10_BEQ     = 4'd10;

    parameter LW    = 7'b0000011,
              SW    = 7'b0100011,
              RTYPE = 7'b0110011,
              ITYPE = 7'b0010011,
              JAL   = 7'b1101111,
              BEQ   = 7'b1100011;

    reg [3:0] current_state, next_state;

    always @(posedge clk or posedge rst) begin
        if (rst)
            current_state <= S0_FETCH;
        else
            current_state <= next_state;
    end

    always @(*) begin
        // Default control values
        Branch     = 0;
        PCUpdate   = 0;
        PCWrite    = 0;
        RegWrite   = 0;
        MemWrite   = 0;
        IRWrite    = 0;
        AdrSrc     = 0;
        ResultSrc  = 2'b00;
        ALUSrcA    = 2'b00;
        ALUSrcB    = 2'b00;
        ALUOp      = 2'b00;
        next_state = current_state;

        case (current_state)
            S0_FETCH: begin
                AdrSrc    = 0;
                IRWrite   = 1;
                ALUSrcA   = 2'b00;
                ALUSrcB   = 2'b10;
                ALUOp     = 2'b00;
                ResultSrc = 2'b10;
                PCUpdate  = 1;
                next_state = S1_DECODE;
            end

            S1_DECODE: begin
                ALUSrcA = 2'b01;
                ALUSrcB = 2'b01;
                ALUOp   = 2'b00;
                case (op)
                    LW:    next_state = S2_MEMADR;
                    SW:    next_state = S2_MEMADR;
                    RTYPE: next_state = S6_EXECUTER;
                    ITYPE: next_state = S8_EXECUTEI;
                    JAL:   next_state = S9_JAL;
                    BEQ:   next_state = S10_BEQ;
                    default: next_state = S0_FETCH;
                endcase
            end

            S2_MEMADR: begin
                ALUSrcA = 2'b10;
                ALUSrcB = 2'b01;
                ALUOp   = 2'b00;
                if (op == LW)
                    next_state = S3_MEMREAD;
                else if (op == SW)
                    next_state = S5_MEMWRITE;
                else
                    next_state = S0_FETCH;
            end

            S3_MEMREAD: begin
                AdrSrc    = 1;
                ResultSrc = 2'b00;
                next_state = S4_MEMWB;
            end

            S4_MEMWB: begin
                ResultSrc = 2'b01;
                RegWrite  = 1;
                next_state = S0_FETCH;
            end

            S5_MEMWRITE: begin
                AdrSrc   = 1;
                MemWrite = 1;
                next_state = S0_FETCH;
            end

            S6_EXECUTER: begin
                ALUSrcA = 2'b10;
                ALUSrcB = 2'b00;
                ALUOp   = 2'b10;
                next_state = S7_ALUWB;
            end

            S7_ALUWB: begin
                ResultSrc = 2'b00;
                RegWrite  = 1;
                next_state = S0_FETCH;
            end

            S8_EXECUTEI: begin
                ALUSrcA = 2'b10;
                ALUSrcB = 2'b01;
                ALUOp   = 2'b10;
                next_state = S7_ALUWB;
            end

            S9_JAL: begin
                ALUSrcA   = 2'b01;
                ALUSrcB   = 2'b10;
                ALUOp     = 2'b00;
                ResultSrc = 2'b00;
                PCUpdate  = 1;
                next_state = S7_ALUWB;
            end

            S10_BEQ: begin
                ALUSrcA   = 2'b10;
                ALUSrcB   = 2'b00;
                ALUOp     = 2'b01;
                Branch    = 1;
                next_state = S0_FETCH;
            end
            
            default: next_state = S0_FETCH;
        endcase
    end

    // PCWrite logic
    always @(*) begin
        PCWrite = (zero & Branch) | PCUpdate;
    end
endmodule

// Control Unit Top
module control_unittop (
    input clk,
    input rst,
    input [6:0] op,
    input [2:0] funct3,
    input funct7,funct7_1,
    input zero,
    output wire PCWrite,
    output wire AdrSrc,
    output wire MemWrite,
    output wire IRWrite,
    output wire [1:0] ResultSrc,
    output wire [3:0] ALUControl,
    output wire [1:0] ALUSrcA,
    output wire [1:0] ALUSrcB,
    output wire [1:0] ImmSrc,
    output wire RegWrite
);

    wire Branch, PCUpdate;
    wire [1:0] ALUOp;
    wire opb5 = op[5];
    wire funct7b5 = funct7;
    wire funct7b1 = funct7_1;

    // FSM instantiation
    main_fsm fsm (
        .clk(clk),
        .rst(rst),
        .zero(zero),
        .op(op),
        .Branch(Branch),
        .PCUpdate(PCUpdate),
        .PCWrite(PCWrite),
        .RegWrite(RegWrite),
        .MemWrite(MemWrite),
        .IRWrite(IRWrite),
        .AdrSrc(AdrSrc),
        .ResultSrc(ResultSrc),
        .ALUSrcA(ALUSrcA),
        .ALUSrcB(ALUSrcB),
        .ALUOp(ALUOp)
    );

    // Immediate Decoder
    Instr_Decoder decoder (
        .op(op),
        .ImmSrc(ImmSrc)
    );

    // ALU Decoder
    aludec alu_decoder (
        .ALUOp(ALUOp),
        .funct3(funct3),
        .opb5(opb5),
        .funct7b5(funct7b5),.funct7b1(funct7b1),
        .ALUControl(ALUControl)
    );
endmodule

// RV32I Multicycle Top Module
module RV32IM_MULTICYCLE(
    input clk,
    input rst,
    output [31:0] Instr,
    output [31:0] SrcA,
    output [31:0] SrcB,
    output [31:0] ALUOut,
    output [31:0] Result,
    output [1:0] ResultSrc,
    output [1:0] ALUSrcA,
    output [1:0] ALUSrcB,
    output [1:0] ImmSrc,
    output [31:0] RD1,
    output [31:0] RD2,
    output [4:0] rs1,
    output [4:0] rs2,
    output [4:0] rd
);

    // Internal wires
    wire PCWrite, AdrSrc, MemWrite, IRWrite, RegWrite;
    wire [3:0] ALUControl;
    wire Zero;
    wire [31:0] PC, OldPC, Adr;
    wire [31:0] RD, A, WriteData, Data, ALUResult, ImmExt;

    // Instruction field extraction
    assign rs1 = Instr[19:15];
    assign rs2 = Instr[24:20];
    assign rd  = Instr[11:7];

    // Module instantiations
    Program_Counter pc(
        .clk(clk),
        .rst(rst),
        .PCWrite(PCWrite),
        .PC_next(Result),
        .PC(PC)
    );

    mux2_1 adr_mux(
        .PC(PC),
        .Result(Result),
        .AdrSrc(AdrSrc),
        .Adr(Adr)
    );

    Instruction_Memory mem(
        .clk(clk),
        .rst(rst),
        .MemWrite(MemWrite),
        .A(Adr),
        .WD(WriteData),
        .RD(RD)
    );

    instruction_register ir(
        .clk(clk),
        .rst(rst),
        .PC(PC),
        .IRWrite(IRWrite),
        .RD(RD),
        .Instr(Instr),
        .OldPC(OldPC)
    );

    Data_Register data_reg(
        .clk(clk),
        .rst(rst),
        .ReadData(RD),
        .Data(Data)
    );

    Register_file rf(
        .clk(clk),
        .rst(rst),
        .we3(RegWrite),
        .A1(rs1),
        .A2(rs2),
        .A3(rd),
        .RD1(RD1),
        .RD2(RD2),
        .WD3(Result)
    );

    rd1_register rd1(
        .clk(clk),
        .rst(rst),
        .RD1(RD1),
        .RD2(RD2),
        .A(A),
        .WriteData(WriteData)
    );

    extend ext(
        .instr(Instr),
        .immsrc(ImmSrc),
        .immext(ImmExt)
    );

    mux1 mux_srcA(
        .PC(PC),
        .OldPC(OldPC),
        .A(A),
        .ALUSrcA(ALUSrcA),
        .SrcA(SrcA)
    );

    mux2 mux_srcB(
        .WriteData(WriteData),
        .ImmExt(ImmExt),
        .ALUSrcB(ALUSrcB),
        .SrcB(SrcB)
    );

    ALU alu(
        .srcA(SrcA),
        .srcB(SrcB),
        .ALUControl(ALUControl),
        .ALUResult(ALUResult),
        .Zero(Zero)
    );

    ALU_register alu_reg(
        .clk(clk),
        .rst(rst),
        .ALUResult(ALUResult),
        .ALUOut(ALUOut)
    );

    mux3 mux_result(
        .ALUOut(ALUOut),
        .Data(Data),
        .ALUResult(ALUResult),
        .ResultSrc(ResultSrc),
        .Result(Result)
    );

    control_unittop ctrl(
        .clk(clk),
        .rst(rst),
        .op(Instr[6:0]),
        .funct3(Instr[14:12]),
        .funct7(Instr[30]),.funct7_1(Instr[25]),
        .zero(Zero),
        .PCWrite(PCWrite),
        .AdrSrc(AdrSrc),
        .MemWrite(MemWrite),
        .IRWrite(IRWrite),
        .ResultSrc(ResultSrc),
        .ALUControl(ALUControl),
        .ALUSrcA(ALUSrcA),
        .ALUSrcB(ALUSrcB),
        .ImmSrc(ImmSrc),
        .RegWrite(RegWrite)
    );
endmodule
