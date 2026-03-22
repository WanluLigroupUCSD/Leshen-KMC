---
name: Push memory with code
description: Every git push to Leshen-KMC must include .claude-memory/ and .autopilot/ directories
type: feedback
---

每次 push 到 Leshen-KMC 仓库时，必须同步推送 `.claude-memory/` 和 `.autopilot/` 目录。

**Why:** 用户要求将项目记忆和 autopilot 状态作为仓库的一部分进行版本管理，确保跨会话的上下文不丢失。

**How to apply:** 在 commit 前，将 `/home/reny0b/.claude/projects/-ibex-user-reny0b-zls-KMC/memory/` 下的文件同步到 `Leshen-KMC/.claude-memory/`，同时确保 `.autopilot/` 状态文件已更新。提交时在 commit message 中注明包含记忆更新。
